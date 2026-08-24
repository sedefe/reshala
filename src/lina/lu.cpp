#include <numeric>

#include "reshala/lina/core/operators.h"
#include "reshala/lina/lina.h"

namespace reshala {

LinaResult Lina::Refactor() {
    stats.n_lus++;

    row_perm.resize(m);
    row_perm_inv.resize(m);
    for (Index i = 0; i < m; ++i) {
        row_perm[i] = row_perm_inv[i] = i;
    }
    Lr.Clear();
    Ur.Clear();
    etas.clear();

    for (Index ic = 0; ic < m; ic++) {
        Index ib = basis_->Basis()[ic];
        for (SvIterator el(Ac_.GetCol(ib)); el; ++el) {
            Ur.GetRow(el.index()).Push(ic, el.value());
        }
    }
    stats.total_nnz_b += Ur.GetNnz();

    row_front.resize(m);
    for (Index i = 0; i < m; ++i) row_front[i].clear();
    for (Index i = 0; i < m; ++i) {
        if (Ur.GetRow(i).Empty()) {
            // std::cerr << "Empty row " << i << "\n";
            return LinaResult::kDegenerate;
        }
        row_front[Ur.GetRow(i).indices()[0]].push_back(i);
    }

    for (Index k = 0; k < m; ++k) {
        // Partial pivoting
        Index pivot_row = -1;
        Scalar pivot_val = 0;

        for (auto& j : row_front[k]) j = row_perm_inv[j];
        std::sort(row_front[k].begin(), row_front[k].end());  // Todo avoid

        bool k_is_in_front = false;
        for (Index i : row_front[k]) {
            k_is_in_front |= (i == k);

            const auto& row_i = Ur.GetRow(i);
            assert(row_i.indices()[0] == k && "Row front is broken");

            Scalar val = row_i.values()[0];
            if (std::abs(val) > std::abs(pivot_val)) {
                pivot_val = val;
                pivot_row = i;
            }
        }
        if (IsZero(pivot_val)) {
            // std::cerr << "Small max_abs @col" << k << ": " << pivot_val << "\n";
            return LinaResult::kDegenerate;
        }

        // Swap rows in U, L, and the permutation vector
        if (pivot_row != k) {
            std::swap(Ur.GetRow(k), Ur.GetRow(pivot_row));
            std::swap(Lr.GetRow(k), Lr.GetRow(pivot_row));

            Index k_perm = row_perm[k];
            Index pivot_perm = row_perm[pivot_row];
            std::swap(row_perm[k], row_perm[pivot_row]);
            std::swap(row_perm_inv[k_perm], row_perm_inv[pivot_perm]);
        }
        u_diag[k] = pivot_val;
        Ur.GetRow(k).EraseOffset(0);

        // Eliminate rows below k
        const auto& row_k = Ur.GetRow(k);
        // Индекс строки из фронта, соответствующий выбранной строке. Без пол-литры не разберёшься:
        // - Если k была во фронте, то вне зависимости от того, был ли свап, пропускать будем строку номер k.
        // - Если не была, то выбрана строка номер pivot_row, и пропускаем её.
        Index skip_row = k_is_in_front ? k : pivot_row;
        for (Index i : row_front[k]) {
            if (i == skip_row) continue;

            auto& row_i = Ur.GetRow(i);
            assert(row_i.indices()[0] == k && "Row front is broken");

            Scalar factor = row_i.values()[0] / pivot_val;

            Lr.GetRow(i).Push(k, factor);  // Store multiplier in L

            // row_i = row_i - factor * row_k
            row_i = axpy(-factor, row_k, row_i);
            row_i.EraseOffset(0);
            if (row_i.Empty()) {
                // std::cerr << "Empty row " << i << "\n";
                return LinaResult::kDegenerate;
            }
            row_front[row_i.indices()[0]].push_back(row_perm[i]);
        }
    }

    Srm2Scm(Lr, Lc);
    Srm2Scm(Ur, Uc);

    stats.total_nnz_l += Lr.GetNnz();
    stats.total_nnz_u += Ur.GetNnz() + m;

    return LinaResult::kOk;
}

}  // namespace reshala
