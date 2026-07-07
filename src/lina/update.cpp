#include "reshala/lina/lina.h"

namespace reshala {

void Lina::Update(Index iv_leaving, Index iv_entering) {
    // Тут мы считаем, что базис уже обновлён
    n_updates_++;

    switch (ut) {
        case UpdType::kSlu:
            Refactor();
            stats.n_lus++;
            break;
        case UpdType::kSluPf:
            if (n_updates_ % kMaxUpdates == 0) {
                Refactor();
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
