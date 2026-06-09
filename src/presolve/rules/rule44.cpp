#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule44::Apply(ModelTracker& tracker) {
    const MilpModel& model = tracker.GetModel();
    Index n_reduced = 0;

    // Todo: count locks?
    for (Index iv = 0; iv < model.GetNVars(); iv++) {
        if (tracker.GetVarMask(iv)) continue;

        const Bounds& bnd = model.GetBounds(iv);
        bool eligible_up = model.GetObj().coefficients[iv] <= 0 and bnd.ri != kInf;
        bool eligible_down = model.GetObj().coefficients[iv] >= 0 and bnd.le != -kInf;

        for (SvIterator el(model.GetCol(iv)); el & (eligible_up || eligible_down); ++el) {
            if (tracker.GetConMask(el.index())) continue;
            const Bounds& rhs = model.GetRhs(el.index());
            if (rhs.le == -kInf) {
                eligible_down &= el.value() >= 0;
                eligible_up &= el.value() <= 0;
            }
            if (rhs.ri == kInf) {
                eligible_down &= el.value() <= 0;
                eligible_up &= el.value() >= 0;
            }
            if (rhs.le != -kInf and rhs.ri != kInf) {
                eligible_down = eligible_up = false;
            }
        }

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
