#include "reshala/lina/lina.h"

#include "reshala/lina/core/operators.h"

namespace reshala {

void Lina::InitBinv() {
    Binv_.ResizeAsZero(m, m);
    for (Index iv = 0; iv < m; iv++) {
        Binv_[iv][iv] = 1;
    }
    n_updates_ = 0;
}

void Lina::Btran(Index iv, DenseVector& res) { res.assign(Binv_[iv], Binv_[iv] + m); }

void Lina::Ftran(Index iv, DenseVector& res) {
    if (iv < n) {
        MulDmSv(Binv_, Ac_->GetCol(iv), res);
    } else {
        MulDmSv(Binv_, SparseVector(m, iv - n, 1.0), res);
    }
}

void Lina::Update(Index iv_leaving, Index iv_entering) {
    // Тут мы считаем, что базис уже обновлён
    n_updates_++;

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

bool Lina::Reinvert() {
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

}  // namespace reshala
