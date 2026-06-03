#include "reshala/lina/lina.h"

#include "reshala/lina/operators.h"

namespace reshala {

void Lina::Init() {
    Binv_.ResizeAsZero(m_, m_);
    for (Index iv = 0; iv < m_; iv++) {
        Binv_.RowView(iv)[iv] = 1;
    }
}

void Lina::Btran(Index iv, DenseVector& res) {
    res.assign(Binv_.RowView(iv), Binv_.RowView(iv) + m_);
}

void Lina::Ftran(const SparseVector& sv, DenseVector& res) { MulDmSv(Binv_, sv, res); }

void Lina::Update(Index iv, const SparseVector& delta) {
    const DenseVector row(Binv_.RowView(iv), Binv_.RowView(iv) + m_);
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
