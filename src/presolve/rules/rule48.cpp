#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule48::Apply(ModelTracker& tracker) {
    const MilpModel& model = tracker.GetModel();
    Index n_reduced = 0;

    for (Index iv = 0; iv < model.GetNVars(); iv++) {
        if (tracker.GetVarMask(iv)) continue;
        if (!model.GetIntegrality(iv)) continue;
        const Bounds& bnd = model.GetBounds(iv);
        if (bnd.le == 0.0) continue;
        if (bnd.ri - bnd.le != 1.0) continue;

        tracker.ConstShiftVar(iv, -bnd.le);
        n_reduced++;
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
