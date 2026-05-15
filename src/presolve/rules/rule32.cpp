#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule32::Apply(ModelInfo& info, std::vector<std::unique_ptr<Transform>>& transforms) {
    MilpModel& model = info.GetModel();
    Index n_reduced = 0;

    for (Index ic = 0; ic < model.GetNCons(); ic++) {
        if (info.GetConMask(ic)) continue;

        const Bounds& act = info.GetActivity(ic);
        const Bounds& rhs = model.GetRhs(ic);
        for (SvIterator el(model.GetAr().GetRow(ic)); el; ++el) {
            if (IsZero(el.value())) continue;

            const Bounds& bnd = model.GetBounds(el.index());
            const Bounds act1 =
                el.value() > 0 ? Bounds{act.le - el.value() * bnd.le, act.ri - el.value() * bnd.ri}
                               : Bounds{act.le - el.value() * bnd.ri, act.ri - el.value() * bnd.le};

            Scalar le_derived, ri_derived;
            if (el.value() > 0) {
                le_derived = (rhs.le - act1.ri) / el.value();
                ri_derived = (rhs.ri - act1.le) / el.value();
            } else {
                le_derived = (rhs.ri - act1.le) / el.value();
                ri_derived = (rhs.le - act1.ri) / el.value();
            }
            if (model.GetIntegrality(el.index())) {
                le_derived = Ceil(le_derived);
                ri_derived = Floor(ri_derived);
            }

            if (StrongGt(le_derived, bnd.le) or StrongLt(ri_derived, bnd.ri)) {
                Bounds new_bnd = {std::max(bnd.le, le_derived), std::min(bnd.ri, ri_derived)};
                info.UpdVarBounds(el.index(), new_bnd);
                n_reduced++;
            }
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
