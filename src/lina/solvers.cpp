#include <cstring>

#include "reshala/lina/core/operators.h"
#include "reshala/lina/lina.h"

namespace reshala {

void Lina::Btran(Index iv, DenseVector& res) {
    DenseVector res_d(m), res_s(m);

    BtranD(iv, res_d);
    // BtranS(iv, res_s);

    res = res_d;

    if (0) {  // Compare
        Scalar d = 0;
        printf("After %d updates:\n", n_updates_);
        for (Index i = 0; i < m; i++) {
            if (!IsZero(res_d[i]) or !IsZero(res_s[i])) {
                printf("[%2d] d %6.2f s %6.2f | d=%5.2f\n", i, res_d[i], res_s[i],
                       res_d[i] - res_s[i]);
            }
            d += std::abs(res_d[i] - res_s[i]);
        }
        printf("Sum d: %.5g\n", d);
    }
}

void Lina::BtranD(Index iv, DenseVector& res) { res.assign(Binv_[iv], Binv_[iv] + m); }

void Lina::BtranS(Index iv, DenseVector& res) {
    // x^T P^T L U = e^T => e = U^T L^T P x

    DenseVector y(m);
    SolveUt(iv, y);  // e = U^T y

    DenseVector x;
    SolveLt(y, x);  // y = L^T x

    for (Index i = 0; i < m; ++i) {  // permute x
        res[perm[i]] = x[i];
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

    // То же, только работаем с помощью разреженной строки
    // SparseVector r(m);
    // const auto& row = Ur.GetRow(iv);
    // for (Index i = 1; i < row.Size(); i++) {
    //     r.Push(row.indices()[i], row.values()[i] * y[iv]);
    // }

    // while (r.Size()) {
    //     Index k = r.indices()[0];
    //     Scalar u_ik = r.values()[0];
    //     Scalar u_kk = Ur.GetRow(k).At(k);

    //     y[k] = -u_ik / u_kk;
    //     r = r - (u_ik / u_kk) * Ur.GetRow(k);
    // }
}

void Lina::SolveLt(const DenseVector& y, DenseVector& x) {
    x = y;
    for (Index ir = m - 1; ir >= 0; ir--) {
        for (SvIterator el(Lc.GetCol(ir)); el; ++el) {
            x[ir] -= el.value() * x[el.index()];
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
