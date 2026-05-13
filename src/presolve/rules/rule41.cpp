#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule41::Apply(ModelInfo& info, std::vector<std::unique_ptr<Transform>>& transforms) {
    MilpModel& model = info.GetModel();
    Index n_reduced = 0;

    // Обратный порядок тут нужен, чтобы в постсолве индексация не поехала
    for (Index iv = model.GetNVars() - 1; iv >= 0; iv--) {
        if (info.GetVarMask(iv)) continue;

        if (model.GetType(iv) == VarType::kFixed) {
            const Bounds& bnd = model.GetBounds(iv);
            Scalar value = (bnd.le + bnd.ri) / 2;
            transforms.push_back(
                std::make_unique<FixVariableTransform>(FixVariableTransform(iv, value)));
            info.FixVar(iv, value);
            info.MaskVar(iv);

            n_reduced++;
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
