#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule41::Apply(ModelInfo& info, std::vector<std::unique_ptr<Transform>>& transforms) {
    MilpModel& model = info.GetModel();
    Index n_reduced = 0;

    for (Index iv = 0; iv < model.GetNVars(); iv++) {
        if (info.GetVarMask(iv)) continue;

        if (model.GetType(iv) == BndType::kFixed) {
            const Bounds& bnd = model.GetBounds(iv);
            Scalar value = (bnd.le + bnd.ri) / 2;

            transforms.push_back(std::make_unique<FixVariableTransform>(
                FixVariableTransform(info.GetOrigVarIdx()[iv], value)));
            info.FixVar(iv, value);
            info.MaskVar(iv);

            n_reduced++;
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
