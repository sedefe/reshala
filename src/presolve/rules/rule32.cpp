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

            Index iv = el.index();
            Scalar val = el.value();
            Bounds bnd = model.GetBounds(iv);
            Bounds derived = tracker.DeriveBounds(ic, iv, tracker.GetActivity(ic), bnd, val);

            if (StrongGt(derived.le, bnd.le) or StrongLt(derived.ri, bnd.ri)) {
                Bounds new_bnd = {std::max(bnd.le, derived.le), std::min(bnd.ri, derived.ri)};
                if (StrongGt(new_bnd.le, new_bnd.ri)) return RuleResult::kInfeasible;

                if (GoodBndChange(bnd, new_bnd)) {
                    tracker.UpdVarBounds(iv, std::move(new_bnd));
                    n_reduced++;
                }
            }
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
