#include "reshala/lina/lina.h"

namespace reshala {

void Lina::Update(Index iv_leaving, Index iv_entering) {
    // Тут мы считаем, что базис уже обновлён
    switch (ut) {
        case UpdType::kSlu:
            Refactor();
            break;
        case UpdType::kSluPf:
            ProdForm(iv_leaving, iv_entering);
            break;
        default:
            break;
    }
}

void Lina::ProdForm(Index iv_leaving, Index iv_entering) {
    stats.n_updates++;
    etas.push_back(Eta(ftran_res, iv_leaving));
}

}  // namespace reshala
