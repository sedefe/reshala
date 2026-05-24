#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule33::Apply(ModelTracker& tracker,
                         std::vector<std::unique_ptr<Transform>>& transforms) {
    const MilpModel& model = tracker.GetModel();
    Index n_reduced = 0;

    for (Index ic = 0; ic < model.GetNCons(); ic++) {
        if (tracker.GetConMask(ic)) continue;
        const Bounds& rhs = model.GetRhs(ic);
        // if (rhs.le != -kInf && rhs.ri != kInf) continue;
        if (std::isinf(rhs.le) + std::isinf(rhs.ri) != 1) continue;
        Bounds& act = tracker.GetActivity(ic);

        for (SvIterator el(model.GetAr().GetRow(ic)); el; ++el) {
            Index iv = el.index();
            if (!model.GetIntegrality(iv)) continue;
            if (tracker.GetVarMask(iv)) continue;
            const Bounds& bnd = model.GetBounds(iv);
            if (bnd.le != 0.0 or bnd.ri != 1) continue;  // Todo handle general integers
            Scalar val = el.value();
            const Bounds act1 = (val >= 0) ? Bounds{act.le - val * bnd.le, act.ri - val * bnd.ri}
                                           : Bounds{act.le - val * bnd.ri, act.ri - val * bnd.le};
            Scalar d = 0;
            if (rhs.le == -kInf) {  // Todo вынести из цикла
                if (val > 0) {
                    d = (rhs.ri - act1.ri);
                    if (val >= d and StrongGt(d, 0)) {
                        tracker.UpdCoeff(ic, iv, val - d);
                        tracker.UpdRhs(ic, {-kInf, act1.ri});
                        n_reduced++;
                    }
                } else {
                    d = (rhs.ri - val - act1.ri);
                    if (-val >= d and StrongGt(d, 0)) {
                        tracker.UpdCoeff(ic, iv, val + d);
                        n_reduced++;
                    }
                }
            } else {
                if (val > 0) {
                    d = (act1.le + val - rhs.le);
                    if (val >= d and StrongGt(d, 0)) {
                        tracker.UpdCoeff(ic, iv, val - d);
                        n_reduced++;
                    }
                } else {
                    d = (act1.le - rhs.le);
                    if (-val >= d and StrongGt(d, 0)) {
                        tracker.UpdCoeff(ic, iv, val + d);
                        tracker.UpdRhs(ic, {act1.le, kInf});
                        n_reduced++;
                    }
                }
            }
        }
    }
    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
