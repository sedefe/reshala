#include <cstring>

#include "reshala/lina/core/operators.h"
#include "reshala/lina/lina.h"

namespace reshala {

void Lina::Btran(Index iv, DenseVector& res) {
    switch (ut) {
        case UpdType::kDluSm:
        case UpdType::kSluSm:
            BtranD(iv, res);
            break;
        case UpdType::kSlu:
            BtranS(iv, res);
            break;
        default:
            break;
    }
}

void Lina::Ftran(Index iv, DenseVector& res) {
    switch (ut) {
        case UpdType::kDluSm:
        case UpdType::kSluSm: {
            FtranD(iv, res);
        } break;
        case UpdType::kSlu:
            FtranS(iv, res);
            break;
        default:
            break;
    }
}

void Lina::BtranD(Index iv, DenseVector& res) { res.assign(Binv_[iv], Binv_[iv] + m); }

void Lina::BtranS(Index iv, DenseVector& res) {
    // x^T P^T L U = e^T => e = U^T L^T P x

    DenseVector x(m);
    x[iv] = 1.0;

    SolveUt(x);  // e = U^T y
    SolveLt(x);  // y = L^T x

    Index i = 0;
    for (Index i = 0; i < m; ++i) {  // permute x
        res[i] = x[row_perm_inv[i]];
    }
}

void Lina::SolveUt(DenseVector& x) {
    for (Index k = 0; k < m; k++) {
        if (IsZero(x[k])) continue;
        Scalar u_ik = x[k];
        Scalar u_kk = Ur.GetRow(k).At(k);
        Scalar factor = u_ik / u_kk;

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

void Lina::FtranD(Index iv, DenseVector& res) {
    if (iv < n) {
        MulDmSv(Binv_, Ac_->GetCol(iv), res);
    } else {
        MulDmSv(Binv_, SparseVector(m, iv - n, 1.0), res);
    }
}

void Lina::FtranS(Index iv, DenseVector& res) {
    // P^T L U x = b => L U x = P b

    res.assign(m, 0.0);
    if (iv < n) {
        for (SvIterator el(Ac_->GetCol(iv)); el; ++el) {
            res[row_perm_inv[el.index()]] = el.value();
        }
    } else {
        res[row_perm_inv[iv - n]] = 1.0;
    }

    SolveL(res);  // b = L y
    SolveU(res);  // y = U x
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
        Scalar u_kk = Ur.GetRow(k).At(k);
        Scalar factor = u_ik / u_kk;

        for (SvIterator el(Uc.GetCol(k)); el; ++el) {
            x[el.index()] -= factor * el.value();
        }
        x[k] = factor;
    }
}

}  // namespace reshala
