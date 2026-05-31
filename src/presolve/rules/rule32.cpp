#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule32::Apply(ModelTracker& tracker) {
    const MilpModel& model = tracker.GetModel();
    Index n_reduced = 0;

    for (Index ic = 0; ic < model.GetNCons(); ic++) {
        if (tracker.GetConMask(ic)) continue;

        const Activity& act = tracker.GetActivity(ic);
        const Bounds& rhs = model.GetRhs(ic);
        for (SvIterator el(model.GetRow(ic)); el; ++el) {
            if (IsZero(el.value())) continue;

            const Bounds& bnd = model.GetBounds(el.index());
            Scalar val = el.value();
            Activity act1 = act;
            act1.RmTerm(val, bnd);
            const Bounds& lhs = act1.GetRange();

            Scalar le_derived, ri_derived;
            if (el.value() > 0) {
                le_derived = (rhs.le - lhs.ri) / val;
                ri_derived = (rhs.ri - lhs.le) / val;
            } else {
                le_derived = (rhs.ri - lhs.le) / val;
                ri_derived = (rhs.le - lhs.ri) / val;
            }
            if (model.GetIntegrality(el.index())) {
                le_derived = Ceil(le_derived);
                ri_derived = Floor(ri_derived);
            }

            if (StrongGt(le_derived, bnd.le) or StrongLt(ri_derived, bnd.ri)) {
                Bounds new_bnd = {std::max(bnd.le, le_derived), std::min(bnd.ri, ri_derived)};
                tracker.UpdVarBounds(el.index(), new_bnd);
                n_reduced++;

                if (le_derived > ri_derived) return RuleResult::kInfeasible;
            }
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
