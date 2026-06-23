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

    DenseVector b(m);
    b[iv] = 1.0;

    DenseVector y(m);
    SolveUt(b, y);  // e = U^T y

    DenseVector x(m);
    SolveLt(y, x);  // y = L^T x

    Index i = 0;
    for (Index i = 0; i < m; ++i) {  // permute x
        res[i] = x[row_perm_inv[i]];
    }
}

void Lina::SolveUt(const DenseVector& b, DenseVector& y) {
    y = b;
    for (Index k = 0; k < m; k++) {
        if (IsZero(y[k])) continue;
        Scalar u_ik = y[k];
        Scalar u_kk = Ur.GetRow(k).At(k);
        Scalar factor = u_ik / u_kk;

        for (SvIterator el(Ur.GetRow(k)); el; ++el) {
            y[el.index()] -= factor * el.value();
        }
        y[k] = factor;
    }
}

void Lina::SolveLt(const DenseVector& y, DenseVector& x) {
    x = y;
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

    DenseVector b(m);
    if (iv < n) {
        for (SvIterator el(Ac_->GetCol(iv)); el; ++el) {
            b[row_perm_inv[el.index()]] = el.value();
        }
    } else {
        b[row_perm_inv[iv - n]] = 1.0;
    }

    DenseVector y(m);
    SolveL(b, y);  // b = L y

    res.assign(m, 0);
    SolveU(y, res);  // y = U x
}

void Lina::SolveL(const DenseVector& b, DenseVector& y) {
    y = b;
    for (Index k = 0; k < m; k++) {
        if (IsZero(y[k])) continue;
        for (SvIterator el(Lc.GetCol(k)); el; ++el) {
            y[el.index()] -= el.value() * y[k];
        }
    }
}

void Lina::SolveU(const DenseVector& y, DenseVector& x) {
    x = y;
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
