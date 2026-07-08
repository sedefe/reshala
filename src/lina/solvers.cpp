#include <cstring>

#include "reshala/lina/core/operators.h"
#include "reshala/lina/lina.h"

namespace reshala {

void Lina::Btran(Index iv, DenseVector& res) {
    // x^T B = ep^T
    // x^T R Bs Cb = ep^T
    // x^T R P^T L U Cb = ep^T
    // x = R^-1 P^T L^-T U^-T E1^-T ... Ek^-T Cb^-1 ep

    DenseVector x(m);
    x[iv] = std::ldexp(1.0, -scaling.col[basis_->Basis()[iv]]);

    for (Index i = etas.size() - 1; i >= 0; i--) {
        EtaBtran(etas[i], x);
    }

    SolveUt(x);  // e = U^T y
    SolveLt(x);  // y = L^T x

    for (Index i = 0; i < m; ++i) {  // permute x
        res[i] = std::ldexp(x[row_perm_inv[i]], -scaling.row[i]);
    }
}

void Lina::Btran(DenseVector& x, DenseVector& res) {
    // x^T B = b^T
    // x^T R Bs Cb = b^T
    // x^T R P^T L U Cb = b^T
    // x = R^-1 P^T L^-T U^-T E1^-T ... Ek^-T Cb^-1 b

    for (Index i = 0; i < m; i++) {
        x[i] = std::ldexp(x[i], -scaling.col[basis_->Basis()[i]]);
    }

    for (Index i = etas.size() - 1; i >= 0; i--) {
        EtaBtran(etas[i], x);
    }

    SolveUt(x);  // e = U^T y
    SolveLt(x);  // y = L^T x

    for (Index i = 0; i < m; ++i) {  // permute x
        res[i] = std::ldexp(x[row_perm_inv[i]], -scaling.row[i]);
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

void Lina::Ftran(Index iv, DenseVector& res) {
    // B x = A eq
    // R Bs Cb x = R As C eq
    // Bs Cb x = As C eq
    // P^T L U E1 .. Ek Cb x = As C eq
    // x = Cb^-1 Ek^-1 .. E1^-1 U^-1 L^-1 P As C eq

    // Ek  = I + (aq - ep) ap^T
    // Eks = I + Cb (aq - ep) ap^T Cb^T

    // Assign b
    res.assign(m, 0.0);
    for (SvIterator el(Ac_.GetCol(iv)); el; ++el) {
        res[row_perm_inv[el.index()]] = std::ldexp(el.value(), scaling.col[iv]);
    }

    SolveL(res);  // b = L y
    SolveU(res);  // y = U x

    for (Index i = 0; i < etas.size(); i++) {
        EtaFtran(etas[i], res);
    }

    for (Index i = 0; i < m; i++) {  // Todo use sparse ftran_res for this loop
        res[i] = std::ldexp(res[i], -scaling.col[basis_->Basis()[i]]);
    }

    ftran_res = res;  // For update
    for (MutableSvIterator el(ftran_res); el; ++el) {
        Index scale = scaling.col[basis_->Basis()[el.index()]] - scaling.col[iv];
        el.valueRef() = std::ldexp(el.value(), scale);
    }
}

void Lina::Ftran(const DenseVector& x, DenseVector& res) {
    // B x = b
    // R Bs Cb x = b
    // R P^T L U E1 .. Ek Cb x = b
    // x = Cb^-1 Ek^-1 .. E1^-1 U^-1 L^-1 P R^-1 b

    // P As C eq
    // res[row_perm_inv[el.index()]] = std::ldexp(el.value(), scaling.col[iv]);
    // P R^-1 b
    for (Index i = 0; i < m; i++) {
        res[row_perm_inv[i]] = std::ldexp(x[i], -scaling.row[i]);
    }

    SolveL(res);  // b = L y
    SolveU(res);  // y = U x

    for (Index i = 0; i < etas.size(); i++) {
        EtaFtran(etas[i], res);
    }

    for (Index i = 0; i < m; i++) {
        res[i] = std::ldexp(res[i], -scaling.col[basis_->Basis()[i]]);
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
