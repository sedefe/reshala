#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule41::Apply(MilpModel& model, std::vector<std::unique_ptr<Transform>>& transforms_) {
    bool reduced = false;
    for (Index iv = 0; iv < model.GetNVars(); iv++) {
        if (model.GetType(iv) == VarType::kFixed) {
            const Bounds& bnd = model.GetBounds(iv);
            Scalar value = (bnd.le + bnd.ri) / 2;
            model.GetObj().c0 += model.GetObj().coefficients[iv] * value;

            FixVariableTransform tr(iv, value);
            tr.Do(model);
            transforms_.push_back(std::make_unique<FixVariableTransform>(tr));

            reduced = true;
        }
    }
    return reduced ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
