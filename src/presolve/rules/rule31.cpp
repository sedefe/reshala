#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule31::Apply(ModelTracker& tracker) {
    const MilpModel& model = tracker.GetModel();
    Index n_reduced = 0;

    for (Index ic = 0; ic < model.GetNCons(); ic++) {
        if (tracker.GetConMask(ic)) continue;

        const Bounds& lhs = tracker.GetConRange(ic);
        const Bounds& rhs = model.GetRhs(ic);
        if (WeakGe(lhs.le, rhs.le) and WeakLe(lhs.ri, rhs.ri)) {
            tracker.MaskCon(ic);
            n_reduced++;
        }

        if (StrongGt(lhs.le, rhs.ri) or StrongLt(lhs.ri, rhs.le)) {
            return RuleResult::kInfeasible;
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
