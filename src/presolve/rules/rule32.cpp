#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule32::Apply(ModelTracker& tracker,
                         std::vector<std::unique_ptr<Transform>>& transforms) {
    const MilpModel& model = tracker.GetModel();
    Index n_reduced = 0;

    for (Index ic = 0; ic < model.GetNCons(); ic++) {
        if (tracker.GetConMask(ic)) continue;

        const Bounds& act = tracker.GetActivity(ic);
        const Bounds& rhs = model.GetRhs(ic);
        for (SvIterator el(model.GetRow(ic)); el; ++el) {
            if (IsZero(el.value())) continue;

            const Bounds& bnd = model.GetBounds(el.index());
            Scalar val = el.value();
            const Bounds act1 = (val >= 0) ? Bounds{act.le - val * bnd.le, act.ri - val * bnd.ri}
                                           : Bounds{act.le - val * bnd.ri, act.ri - val * bnd.le};

            Scalar le_derived, ri_derived;
            if (el.value() > 0) {
                le_derived = (rhs.le - act1.ri) / val;
                ri_derived = (rhs.ri - act1.le) / val;
            } else {
                le_derived = (rhs.ri - act1.le) / val;
                ri_derived = (rhs.le - act1.ri) / val;
            }
            if (model.GetIntegrality(el.index())) {
                le_derived = Ceil(le_derived);
                ri_derived = Floor(ri_derived);
            }

            if (StrongGt(le_derived, bnd.le) or StrongLt(ri_derived, bnd.ri)) {
                Bounds new_bnd = {std::max(bnd.le, le_derived), std::min(bnd.ri, ri_derived)};
                tracker.UpdVarBounds(el.index(), new_bnd);
                n_reduced++;
            }
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
