#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule41::Apply(ModelTracker& tracker,
                         std::vector<std::unique_ptr<Transform>>& transforms) {
    const MilpModel& model = tracker.GetModel();
    Index n_reduced = 0;

    for (Index iv = 0; iv < model.GetNVars(); iv++) {
        if (tracker.GetVarMask(iv)) continue;

        if (model.GetType(iv) == BndType::kFixed) {
            const Bounds& bnd = model.GetBounds(iv);
            Scalar value = (bnd.le + bnd.ri) / 2;

            transforms.push_back(std::make_unique<FixVariableTransform>(
                FixVariableTransform(tracker.GetOrigVarIdx()[iv], value)));
            tracker.FixVar(iv, value);
            tracker.MaskVar(iv);

            n_reduced++;
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
