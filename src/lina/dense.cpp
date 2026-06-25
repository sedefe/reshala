#include <cstring>

#include "reshala/lina/core/operators.h"
#include "reshala/lina/lina.h"

namespace reshala {

void Lina::BtranD(Index iv, DenseVector& res) { res.assign(Binv_[iv], Binv_[iv] + m); }

void Lina::FtranD(Index iv, DenseVector& res) {
    if (iv < n) {
        MulDmSv(Binv_, Ac_->GetCol(iv), res);
    } else {
        MulDmSv(Binv_, SparseVector(m, iv - n, 1.0), res);
    }
}

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

bool Lina::SparseLu2Binv() {
    DenseVector e(m), res(m);
    for (Index ir = 0; ir < m; ir++) {
        BtranS(ir, res);
        std::memcpy(Binv_[ir], res.data(), m * sizeof(Scalar));
    }
    return true;
}

void Lina::SherMor(Index iv_leaving, Index iv_entering) {
    Index i_b = basis_->Basis()[iv_leaving];       // Переменная, которая только что вошла в базис
    Index i_nb = basis_->NonBasis()[iv_entering];  // Переменная, которая только что вышла
    const SparseVector& entering_col = i_b < n ? Ac_->GetCol(i_b) : SparseVector(m, i_b - n, 1.0);
    const SparseVector& leaving_col = i_nb < n ? Ac_->GetCol(i_nb) : SparseVector(m, i_nb - n, 1.0);
    const SparseVector delta = entering_col - leaving_col;

    const DenseVector row(Binv_[iv_leaving], Binv_[iv_leaving] + m);
    Scalar multiplier;
    dot(row, delta, multiplier);
    multiplier = 1 / (1 + multiplier);
    SparseVector row_sv(row);
    row_sv *= multiplier;

    DenseVector d(m, 0);
    MulDmSv(Binv_, delta, d);

    for (SvIterator el(row_sv); el; ++el) {
        for (Index i = 0; i < m; i++) {
            Binv_[i][el.index()] -= d[i] * el.value();
        }
    }
}

}  // namespace reshala
