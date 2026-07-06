#include <numeric>

#include "reshala/lina/core/operators.h"
#include "reshala/lina/lina.h"

namespace reshala {

bool Lina::SparseLU() {
    row_perm.resize(m);
    row_perm_inv.resize(m);
    for (Index i = 0; i < m; ++i) {
        row_perm[i] = row_perm_inv[i] = i;
    }
    Lr.Clear();
    Ur.Clear();

    for (Index ic = 0; ic < m; ic++) {
        Index ib = basis_->Basis()[ic];
        for (SvIterator el(Ac_.GetCol(ib)); el; ++el) {
            Ur.GetRow(el.index()).Push(ic, el.value());
        }
    }

    for (Index k = 0; k < m; ++k) {
        // Partial pivoting
        Index pivot_row = k;
        Scalar pivot_val = 0;
        for (Index i = k; i < m; ++i) {
            Scalar val = Ur.GetRow(i).At(k);
            if (std::abs(val) > std::abs(pivot_val)) {
                pivot_val = val;
                pivot_row = i;
            }
        }
        if (IsZero(pivot_val)) {
            std::cerr << "Small max_abs @col" << k << ": " << pivot_val << "\n";
            exit(0);
            return false;
        }

        // Swap rows in U, L, and the permutation vector
        if (pivot_row != k) {
            std::swap(Ur.GetRow(k), Ur.GetRow(pivot_row));
            std::swap(Lr.GetRow(k), Lr.GetRow(pivot_row));

            Index row_k = row_perm[k];
            Index row_pivot = row_perm[pivot_row];
            std::swap(row_perm[k], row_perm[pivot_row]);
            std::swap(row_perm_inv[row_k], row_perm_inv[row_pivot]);
        }
        u_diag[k] = pivot_val;
        Ur.GetRow(k).EraseOffset(0);

        // Eliminate rows below k
        for (Index i = k + 1; i < m; ++i) {
            if (Ur.GetRow(i).indices()[0] != k) continue;  // a_ik is already zero
            Scalar factor = Ur.GetRow(i).values()[0] / pivot_val;

            Lr.GetRow(i).Push(k, factor);  // Store multiplier in L

            // row_i = row_i - factor * row_k
            const auto& row_k = Ur.GetRow(k);
            SparseVector scaled = row_k * factor;  // Todo: combine to daxpy
            Ur.GetRow(i) = Ur.GetRow(i) - scaled;
            Ur.GetRow(i).EraseOffset(0);
        }
    }

    Srm2Scm(Lr, Lc);
    Srm2Scm(Ur, Uc);

    return true;
}

}  // namespace reshala
