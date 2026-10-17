#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule44::Apply(ModelTracker& tracker, Index iv) {
    const MilpModel& model = tracker.GetModel();

    const Bounds& bnd = model.GetBounds(iv);
    bool eligible_down = model.GetObj().coefficients[iv] >= 0 and bnd.le != -kInf and
                         tracker.GetLocks().GetNLocks(iv, LockType::kDown) == 0;
    bool eligible_up = model.GetObj().coefficients[iv] <= 0 and bnd.ri != kInf and
                       tracker.GetLocks().GetNLocks(iv, LockType::kUp) == 0;

    if (eligible_down) {
        tracker.FixVar(iv, bnd.le);
        return RuleResult::kReduced;
    }
    if (eligible_up) {
        tracker.FixVar(iv, bnd.ri);
        return RuleResult::kReduced;
    }

    return RuleResult::kUnchanged;
}

}  // namespace reshala
