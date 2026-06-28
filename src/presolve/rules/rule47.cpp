#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule47::Apply(ModelTracker& tracker) {
    const MilpModel& model = tracker.GetModel();
    Index n_reduced = 0;
    Index size = 0;

    for (Index iv = 0; iv < model.GetNVars(); iv++) {
        if (tracker.GetVarMask(iv)) continue;

        Index ic = -1;
        Scalar a = kNan;
        for (SvIterator el(model.GetCol(iv)); el; ++el) {
            if (tracker.GetConMask(el.index())) continue;
            if (ic != -1) {
                ic = -1;
                break;
            }
            ic = el.index();
            a = el.value();
        }
        if (ic < 0) continue;
        const auto& rhs = model.GetRhs(ic);
        if (!IsZero(rhs.ri - rhs.le)) continue;
        if (model.GetIntegrality(iv)) {
            // Todo
        } else {
            tracker.SlackSub(iv, ic, a);
            n_reduced++;
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
