#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule48::Apply(ModelTracker& tracker, Index iv) {
    const MilpModel& model = tracker.GetModel();

    if (!model.GetIntegrality(iv)) return RuleResult::kUnchanged;
    const Bounds& bnd = model.GetBounds(iv);
    if (bnd.le == 0.0) return RuleResult::kUnchanged;
    if (bnd.ri - bnd.le != 1.0) return RuleResult::kUnchanged;

    tracker.ConstShiftVar(iv, -bnd.le);

    return RuleResult::kReduced;
}

}  // namespace reshala
