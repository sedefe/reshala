#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule44::Apply(ModelTracker& tracker) {
    const MilpModel& model = tracker.GetModel();
    Index n_reduced = 0;

    // Todo: count locks?
    for (Index iv = 0; iv < model.GetNVars(); iv++) {
        if (tracker.VarIsDeleted(iv)) continue;

        const Bounds& bnd = model.GetBounds(iv);
        bool eligible_down = model.GetObj().coefficients[iv] >= 0 and bnd.le != -kInf and
                             tracker.GetLocks().GetNLocks(iv, LockType::kDown) == 0;
        bool eligible_up = model.GetObj().coefficients[iv] <= 0 and bnd.ri != kInf and
                           tracker.GetLocks().GetNLocks(iv, LockType::kUp) == 0;

        if (eligible_down) {
            Scalar value = bnd.le;
            tracker.FixVar(iv, value);
            n_reduced++;
        } else if (eligible_up) {
            Scalar value = bnd.ri;
            tracker.FixVar(iv, value);
            n_reduced++;
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
