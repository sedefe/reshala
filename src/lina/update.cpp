#include <cstring>

#include "reshala/lina/core/operators.h"
#include "reshala/lina/lina.h"

namespace reshala {

void Lina::Update(Index iv_leaving, Index iv_entering) {
    // Тут мы считаем, что базис уже обновлён
    n_updates_++;

    switch (ut) {
        case UpdType::kDluSm:
            if (n_updates_ % kMaxUpdates == 0) {
                InvertD();
            } else {
                SherMor(iv_leaving, iv_entering);
            }
            break;
        case UpdType::kSluSm:
            if (n_updates_ % kMaxUpdates == 0) {
                InvertS();
                BuildBinv();
            } else {
                SherMor(iv_leaving, iv_entering);
            }
            break;
        case UpdType::kSlu:
            InvertS();
            break;
        default:
            break;
    }
}

bool Lina::BuildBinv() {
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
