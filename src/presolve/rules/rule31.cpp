#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule31::Apply(ModelTracker& tracker) {
    const MilpModel& model = tracker.GetModel();
    Index n_reduced = 0;

    for (Index ic = 0; ic < model.GetNCons(); ic++) {
        if (tracker.GetConMask(ic)) continue;

        const Bounds& lhs = tracker.GetConRange(ic);
        const Bounds& rhs = model.GetRhs(ic);

        if (StrongGt(lhs.le, rhs.ri) or StrongLt(lhs.ri, rhs.le)) {
            return RuleResult::kInfeasible;
        }

        if (WeakGe(lhs.le, rhs.le) and rhs.le != -kInf) {  // Release lower
            tracker.UpdRhs(ic, {-kInf, rhs.ri});
            n_reduced++;
        }

        if (WeakLe(lhs.ri, rhs.ri) and rhs.ri != kInf) {  // Release upper
            tracker.UpdRhs(ic, {rhs.le, kInf});
            n_reduced++;
        }

        if (rhs.le == -kInf and rhs.ri == kInf) {
            tracker.MaskCon(ic);
            n_reduced++;
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
