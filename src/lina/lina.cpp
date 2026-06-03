#include "reshala/lina/lina.h"

#include "reshala/lina/core/operators.h"

namespace reshala {

void Lina::InitBinv() {
    Binv_.ResizeAsZero(m_, m_);
    for (Index iv = 0; iv < m_; iv++) {
        Binv_.RowView(iv)[iv] = 1;
    }
}

void Lina::Btran(Index iv, DenseVector& res) {
    res.assign(Binv_.RowView(iv), Binv_.RowView(iv) + m_);
}

void Lina::Ftran(Index iv, DenseVector& res) {
    if (iv < n_) {
        MulDmSv(Binv_, Ac_->GetCol(iv), res);
    } else {
        MulDmSv(Binv_, SparseVector(m_, iv - n_, 1.0), res);
    }
}

void Lina::Update(Index iv_leaving, Index iv_entering) {
    Index i_b = basis_->Basis()[iv_leaving];
    Index i_nb = basis_->NonBasis()[iv_entering];
    const SparseVector& leaving_col = i_b < n_ ? Ac_->GetCol(i_b) : SparseVector(m_, i_b - n_, 1.0);
    const SparseVector& entering_col =
        i_nb < n_ ? Ac_->GetCol(i_nb) : SparseVector(m_, i_nb - n_, 1.0);
    const SparseVector delta = entering_col - leaving_col;

    const DenseVector row(Binv_.RowView(iv_leaving), Binv_.RowView(iv_leaving) + m_);
    Scalar multiplier;
    dot(row, delta, multiplier);
    multiplier = 1 / (1 + multiplier);
    SparseVector row_sv(row);
    row_sv *= multiplier;

    DenseVector d(m_, 0);
    MulDmSv(Binv_, delta, d);

    for (Index i = 0; i < m_; i++) {
        for (SvIterator el(row_sv); el; ++el) {
            Binv_.RowView(i)[el.index()] -= d[i] * el.value();
        }
    }
}

}  // namespace reshala
