#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule32::Apply(ModelTracker& tracker, Index ic) {
    const MilpModel& model = tracker.GetModel();
    Index n_reduced = 0;

    const Activity& act = tracker.GetActivity(ic);
    const Bounds& rhs = model.GetRhs(ic);
    for (SvIterator el(model.GetRow(ic)); el; ++el) {
        Index iv = el.index();

        if (tracker.VarIsDeleted(iv)) continue;
        if (IsZero(el.value())) continue;

        Scalar val = el.value();
        Bounds bnd = model.GetBounds(iv);  // todo: &?
        Bounds derived = tracker.DeriveBounds(ic, iv, tracker.GetActivity(ic), bnd, val);

        if (StrongGt(derived.le, bnd.le) or StrongLt(derived.ri, bnd.ri)) {
            derived = BoundsIntersection(bnd, derived);
            if (StrongGt(derived.le, derived.ri)) return RuleResult::kInfeasible;

            if (GoodBndChange(bnd, derived)) {
                tracker.UpdVarBounds(iv, std::move(derived));
                n_reduced++;
            }
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
