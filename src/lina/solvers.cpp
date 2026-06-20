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

void Lina::BtranD(Index iv, DenseVector& res) { res.assign(Binv_[iv], Binv_[iv] + m); }

void Lina::BtranS(Index iv, DenseVector& res) {
    // x^T P^T L U = e^T => e = U^T L^T P x

    DenseVector y(m);
    SolveUt(iv, y);  // e = U^T y

    DenseVector x(m);
    SolveLt(y, x);  // y = L^T x

    Index i = 0;
    for (Index i = 0; i < m; ++i) {  // permute x
        res[i] = x[row_perm_inv[i]];
    }
}

void Lina::SolveUt(Index iv, DenseVector& y) {
    // Мы знаем, что правая часть системы - просто орт. Поэтому сразу используем соответствующую
    // строку, а потом убираем из неё всё, кроме первой единицы, с помощью следующих строк

    Scalar u_ii = Ur.GetRow(iv).At(iv);
    y[iv] = 1 / u_ii;

    const auto& row = Ur.GetRow(iv);
    for (Index i = 1; i < row.Size(); i++) {
        y[row.indices()[i]] = row.values()[i] * y[iv];
    }

    for (Index k = iv + 1; k < m; k++) {
        if (IsZero(y[k])) continue;
        Scalar u_ik = y[k];
        Scalar u_kk = Ur.GetRow(k).At(k);
        Scalar factor = u_ik / u_kk;

        for (SvIterator el(Ur.GetRow(k)); el; ++el) {
            y[el.index()] -= factor * el.value();
        }
        y[k] = -factor;
    }

    // // То же, только работаем с помощью разреженной строки
    // SparseVector residual(m);
    // const auto& row = Ur.GetRow(iv);
    // for (Index i = 1; i < row.Size(); i++) {
    //     residual.Push(row.indices()[i], row.values()[i] * y[iv]);
    // }

    // while (residual.Size()) {
    //     Index k = residual.indices()[0];
    //     Scalar u_ik = residual.values()[0];
    //     Scalar u_kk = Ur.GetRow(k).At(k);

    //     y[k] = -u_ik / u_kk;
    //     residual = residual - (u_ik / u_kk) * Ur.GetRow(k);
    // }
}

void Lina::SolveLt(DenseVector& y, DenseVector& x) {
    for (Index k = m - 1; k >= 0; k--) {
        if (IsZero(y[k])) continue;
        x[k] = y[k];
        for (SvIterator el(Lr.GetRow(k)); el; ++el) {
            y[el.index()] -= el.value() * x[k];
        }
    }
}

void Lina::Ftran(Index iv, DenseVector& res) {
    if (iv < n) {
        MulDmSv(Binv_, Ac_->GetCol(iv), res);
    } else {
        MulDmSv(Binv_, SparseVector(m, iv - n, 1.0), res);
    }
}

}  // namespace reshala
