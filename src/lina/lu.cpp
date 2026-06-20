#include <cstring>
#include <numeric>

#include "reshala/lina/core/operators.h"
#include "reshala/lina/lina.h"

namespace reshala {

bool Lina::InvertD() {
    Binv_.ResizeAsZero(m, m);
    for (Index iv = 0; iv < m; iv++) {
        Index i_b = basis_->Basis()[iv];
        if (i_b < n) {
            const auto& col = Ac_->GetCol(i_b);
            for (SvIterator el(col); el; ++el) {
                Binv_[el.index()][iv] = el.value();
            }
        } else {
            Binv_[i_b - n][iv] = 1.0;
        }
    }

    return Invert(Binv_);
}

bool Lina::InvertS() {
    SparseRowMatrix B(m, m);
    for (Index ic = 0; ic < m; ic++) {
        Index ib = basis_->Basis()[ic];
        if (ib < n) {
            for (SvIterator el(Ac_->GetCol(ib)); el; ++el) {
                B.GetRow(el.index()).Push(ic, el.value());
            }
        } else
            B.GetRow(ib - n).Push(ic, 1);
    }
    bool lu_res = SparseLU(B);
    if (!lu_res) return lu_res;

    DenseVector e(m), res(m);
    for (Index ir = 0; ir < m; ir++) {
        BtranS(ir, res);
        std::memcpy(Binv_[ir], res.data(), m * sizeof(Scalar));
    }
    return true;
}

bool Lina::SparseLU(const SparseRowMatrix& A) {
    assert(A.GetNRows() == m);
    assert(A.GetNCols() == m);

    Init();
    Ur = A;

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

        // --- Eliminate rows below k ---
        for (Index i = k + 1; i < m; ++i) {
            Scalar aik = Ur.GetRow(i).At(k);
            if (IsZero(aik)) continue;  // already zero

            Scalar factor = aik / pivot_val;

            Lr.GetRow(i).Push(k, factor);  // Store multiplier in L

            // row_i = row_i - factor * row_k
            const auto& row_k = Ur.GetRow(k);
            SparseVector scaled = row_k * factor;  // Todo: combine to daxpy
            Ur.GetRow(i) = Ur.GetRow(i) - scaled;
        }
    }

    Srm2Scm(Lr, Lc);
    Srm2Scm(Ur, Uc);

    return true;
}

}  // namespace reshala
