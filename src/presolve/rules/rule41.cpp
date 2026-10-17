#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule41::Apply(ModelTracker& tracker, Index iv) {
    const MilpModel& model = tracker.GetModel();

    if (model.GetType(iv) == BndType::kFixed) {
        const Bounds& bnd = model.GetBounds(iv);
        Scalar value = (bnd.le + bnd.ri) / 2;
        tracker.FixVar(iv, value);
        return RuleResult::kReduced;
    }

    return RuleResult::kUnchanged;
}

}  // namespace reshala
