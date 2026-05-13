#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule32::Apply(ModelInfo& info, std::vector<std::unique_ptr<Transform>>& transforms) {
    MilpModel& model = info.GetModel();
    Index n_reduced = 0;

    for (Index ic = 0; ic < model.GetNCons(); ic++) {
        if (info.GetConMask(ic)) continue;

        const Bounds& act = info.GetActivity(ic);
        const Bounds& rhs = model.GetRhs()[ic];
        for (SvIterator el(model.GetAr().GetRow(ic)); el; ++el) {
            const Bounds& bnd = model.GetBounds(el.index());
            const Bounds act1 =
                el.value() >= 0
                    ? Bounds{act.le - el.value() * bnd.le, act.ri - el.value() * bnd.ri}
                    : Bounds{act.le - el.value() * bnd.ri, act.ri - el.value() * bnd.le};
            Scalar ri_derived = (rhs.ri - act1.le) / el.value();

            if (StrongGt(el.value(), 0.0)) {
                if (StrongLt(ri_derived, bnd.ri)) {
                    info.UpdVarBounds(el.index(), {bnd.le, ri_derived});
                    n_reduced++;
                }
            } else if (StrongLt(el.value(), 0.0)) {
            }
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
