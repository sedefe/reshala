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
                SparseLu2Binv();
            } else {
                SherMor(iv_leaving, iv_entering);
            }
            break;
        case UpdType::kSlu:
            InvertS();
            break;
        case UpdType::kSluPf:
            if (n_updates_ % kMaxUpdates == 0) {
                InvertS();
                etas.clear();
            } else {
                ProdForm(iv_leaving, iv_entering);
            }
        default:
            break;
    }
}

void Lina::ProdForm(Index iv_leaving, Index iv_entering) {
    etas.push_back(Eta(ftran_res, iv_leaving));
}

}  // namespace reshala
