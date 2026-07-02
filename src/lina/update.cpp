#include "reshala/lina/lina.h"

namespace reshala {

void Lina::Update(Index iv_leaving, Index iv_entering) {
    // Тут мы считаем, что базис уже обновлён
    n_updates_++;

    switch (ut) {
        case UpdType::kDluSm:
            if (n_updates_ % kMaxUpdates == 0) {
                InvertD();
                stats.n_lus++;
            } else {
                SherMor(iv_leaving, iv_entering);
                stats.n_updates++;
            }
            break;
        case UpdType::kSluSm:
            if (n_updates_ % kMaxUpdates == 0) {
                SparseLU();
                SparseLu2Binv();
                stats.n_lus++;
            } else {
                SherMor(iv_leaving, iv_entering);
                stats.n_updates++;
            }
            break;
        case UpdType::kSlu:
            SparseLU();
            stats.n_lus++;
            break;
        case UpdType::kSluPf:
            if (n_updates_ % kMaxUpdates == 0) {
                SparseLU();
                etas.clear();
                stats.n_lus++;
            } else {
                ProdForm(iv_leaving, iv_entering);
                stats.n_updates++;
            }
        default:
            break;
    }
}

void Lina::ProdForm(Index iv_leaving, Index iv_entering) {
    etas.push_back(Eta(ftran_res, iv_leaving));
}

}  // namespace reshala
