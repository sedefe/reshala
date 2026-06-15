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
    // x^T P^T L U = b^T => b = U^T L^T P x

    DenseVector b(m);
    b[iv] = 1;

    DenseVector y(m);
    SolveUt(b, y);  // b = U^T y

    DenseVector x;
    SolveLt(y, x);  // y = L^T x

    for (Index i = 0; i < m; ++i) {  // permute x
        res[perm[i]] = x[i];
    }
}

void Lina::SolveUt(const DenseVector& b, DenseVector& y) {
    for (Index ir = 0; ir < m; ir++) {
        y[ir] = b[ir];
        for (SvIterator el(Uc.GetCol(ir)); el; ++el) {
            if (el.index() != ir)
                y[ir] -= el.value() * y[el.index()];
            else {
                assert(!IsZero(el.value()) && "Degenerate U");
                y[ir] /= el.value();
            }
        }
    }
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
