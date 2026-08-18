#include "reshala/lina/lina.h"

namespace reshala {

LinaResult Lina::Update(Index iv_leaving, Index iv_entering) {
    // Тут мы считаем, что базис уже обновлён
    switch (ut) {
        case UpdType::kLu:
            return Refactor();
            break;
        case UpdType::kPf:
            return ProdForm(iv_leaving, iv_entering);
            break;
        default:
            break;
    }
    return LinaResult::kUnknown;
}

LinaResult Lina::ProdForm(Index iv_leaving, Index iv_entering) {
    stats.n_updates++;
    Eta eta(ftran_res, iv_leaving);
    etas.push_back(eta);
    if (IsZero(eta.diag)) {
        return LinaResult::kDegenerate;
    }
    return LinaResult::kOk;
}

}  // namespace reshala
