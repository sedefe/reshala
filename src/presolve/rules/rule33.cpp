#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule33::Apply(ModelTracker& info, std::vector<std::unique_ptr<Transform>>& transforms) {
    const MilpModel& model = info.GetModel();
    Index n_reduced = 0;

    for (Index ic = 0; ic < model.GetNCons(); ic++) {
        if (info.GetConMask(ic)) continue;
        const Bounds& rhs = model.GetRhs(ic);
        // if (rhs.le != -kInf && rhs.ri != kInf) continue;
        if (std::isinf(rhs.le) + std::isinf(rhs.ri) != 1) continue;
        Bounds& act = info.GetActivity(ic);

        for (SvIterator el(model.GetAr().GetRow(ic)); el; ++el) {
            Index iv = el.index();
            if (!model.GetIntegrality(iv)) continue;
            if (info.GetVarMask(iv)) continue;
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
                        info.UpdCoeff(ic, iv, val - d);
                        info.UpdRhs(ic, {-kInf, act1.ri});
                        n_reduced++;
                    }
                } else {
                    d = (rhs.ri - val - act1.ri);
                    if (-val >= d and StrongGt(d, 0)) {
                        info.UpdCoeff(ic, iv, val + d);
                        // info.UpdRhs(ic, {-kInf, rhs.ri});
                        n_reduced++;
                    }
                }
            } else {
                if (val > 0) {
                    d = (act1.le + val - rhs.le);
                    if (val >= d and StrongGt(d, 0)) {
                        info.UpdCoeff(ic, iv, val - d);
                        // info.UpdRhs(ic, {rhs.le, kInf});
                        n_reduced++;
                    }
                } else {
                    d = (act1.le - rhs.le);
                    if (-val >= d and StrongGt(d, 0)) {
                        info.UpdCoeff(ic, iv, val + d);
                        info.UpdRhs(ic, {act1.le, kInf});
                        n_reduced++;
                    }
                }
            }
        }
    }
    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
