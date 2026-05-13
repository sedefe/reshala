#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule31::Apply(ModelInfo& info, std::vector<std::unique_ptr<Transform>>& transforms) {
    MilpModel& model = info.GetModel();
    Index n_reduced = 0;

    // Обратный порядок тут нужен, чтобы в постсолве индексация не поехала
    for (Index ic = 0; ic < model.GetNCons(); ic++) {
        if (info.GetConMask(ic)) continue;

        const Bounds& act = info.GetActivity(ic);
        const Bounds& rhs = model.GetRhs()[ic];
        if (SoftGe(act.le, rhs.le) and SoftLe(act.ri, rhs.ri)) {
            info.MaskCon(ic);
            n_reduced++;
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
