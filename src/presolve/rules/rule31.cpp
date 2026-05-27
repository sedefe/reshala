#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule31::Apply(ModelTracker& tracker) {
    const MilpModel& model = tracker.GetModel();
    Index n_reduced = 0;

    for (Index ic = 0; ic < model.GetNCons(); ic++) {
        if (tracker.GetConMask(ic)) continue;

        const Bounds& act = tracker.GetActivity(ic);
        const Bounds& rhs = model.GetRhs(ic);
        if (WeakGe(act.le, rhs.le) and WeakLe(act.ri, rhs.ri)) {
            tracker.MaskCon(ic);
            n_reduced++;
        }

        if (StrongGt(act.le, rhs.ri) or StrongLt(act.ri, rhs.le)) {
            tracker.ClaimInfeasible();
            break;
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
