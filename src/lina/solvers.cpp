#include <cstring>

#include "reshala/lina/core/operators.h"
#include "reshala/lina/lina.h"

namespace reshala {

void Lina::Btran(const SparseVector& b, DenseVector& res) {
    // x^T B = b^T
    // x^T P^T L U E1 .. Ek = b^T
    // x = P^T L^-T U^-T E1^-T ... Ek^-T b

    DenseVector x = b.ToDense();

    for (Index i = etas.size() - 1; i >= 0; i--) {
        EtaBtran(etas[i], x);
    }

    SolveUt(x);  // e = U^T y
    SolveLt(x);  // y = L^T x

    for (Index i = 0; i < m; ++i) {  // permute x
        res[i] = x[row_perm_inv[i]];
    }
}

void Lina::Btran(DenseVector& x, DenseVector& res) {
    for (Index i = etas.size() - 1; i >= 0; i--) {
        EtaBtran(etas[i], x);
    }

    SolveUt(x);  // e = U^T y
    SolveLt(x);  // y = L^T x

    for (Index i = 0; i < m; ++i) {  // permute x
        res[i] = x[row_perm_inv[i]];
    }
}

void Lina::SolveUt(DenseVector& x) {
    for (Index k = 0; k < m; k++) {
        if (IsZero(x[k])) continue;
        Scalar u_ik = x[k];
        Scalar factor = u_ik / u_diag[k];

        for (SvIterator el(Ur.GetRow(k)); el; ++el) {
            x[el.index()] -= factor * el.value();
        }
        x[k] = factor;
    }
}

void Lina::SolveLt(DenseVector& x) {
    for (Index k = m - 1; k >= 0; k--) {
        if (IsZero(x[k])) continue;
        for (SvIterator el(Lr.GetRow(k)); el; ++el) {
            x[el.index()] -= el.value() * x[k];
        }
    }
}

void Lina::EtaBtran(const Eta& eta, DenseVector& x) {
    Index p = eta.i_col;
    Scalar d;
    dot(x, eta.eta, d);
    x[p] += (x[p] - d) / eta.diag;
}

void Lina::Ftran(const SparseVector& b, DenseVector& res) {
    // B x = b
    // P^T L U E1 .. Ek x = b
    // x = Ek^-1 .. E1^-1 U^-1 L^-1 P b

    // Ek  = I + (aq - ep) ap^T

    // Assign b
    res.assign(m, 0.0);
    for (SvIterator el(b); el; ++el) {
        res[row_perm_inv[el.index()]] = el.value();
    }

    SolveL(res);  // b = L y
    SolveU(res);  // y = U x

    for (Index i = 0; i < etas.size(); i++) {
        EtaFtran(etas[i], res);
    }

    ftran_res = res;  // For update
}

void Lina::Ftran(const DenseVector& x, DenseVector& res) {
    for (Index i = 0; i < m; i++) {
        res[row_perm_inv[i]] = x[i];
    }

    SolveL(res);  // b = L y
    SolveU(res);  // y = U x

    for (Index i = 0; i < etas.size(); i++) {
        EtaFtran(etas[i], res);
    }
}

void Lina::SolveL(DenseVector& x) {
    for (Index k = 0; k < m; k++) {
        if (IsZero(x[k])) continue;
        for (SvIterator el(Lc.GetCol(k)); el; ++el) {
            x[el.index()] -= el.value() * x[k];
        }
    }
}

void Lina::SolveU(DenseVector& x) {
    for (Index k = m - 1; k >= 0; k--) {
        if (IsZero(x[k])) continue;
        Scalar u_ik = x[k];
        Scalar factor = u_ik / u_diag[k];

        for (SvIterator el(Uc.GetCol(k)); el; ++el) {
            x[el.index()] -= factor * el.value();
        }
        x[k] = factor;
    }
}

void Lina::EtaFtran(const Eta& eta, DenseVector& x) {
    Index p = eta.i_col;
    Scalar xp = x[p] / eta.diag;
    for (SvIterator el(eta.eta); el; ++el) {
        x[el.index()] -= xp * el.value();
    }
    x[p] = xp;
}

}  // namespace reshala
