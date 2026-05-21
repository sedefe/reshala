#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule33::Apply(ModelInfo& info, std::vector<std::unique_ptr<Transform>>& transforms) {
    MilpModel& model = info.GetModel();
    Index n_reduced = 0;

    for (Index ic = 0; ic < model.GetNCons(); ic++) {
        if (info.GetConMask(ic)) continue;
        const Bounds& rhs = model.GetRhs(ic);
        if (rhs.le != -kInf && rhs.ri != kInf) continue;
        Bounds& act = info.GetActivity(ic);

        for (SvIterator el(model.GetAr().GetRow(ic)); el; ++el) {
            Index iv = el.index();
            if (!model.GetIntegrality(iv)) continue;
            if (info.GetVarMask(iv)) continue;
            const Bounds& bnd = model.GetBounds(iv);
            if (bnd.le != 0.0 or bnd.ri != 1) continue;
            Scalar val = el.value();
            const Bounds act1 = Bounds{act.le - val * bnd.le, act.ri - val * bnd.ri};

            if (val > 0) {
                if (rhs.le == -kInf) { // Todo вынести из цикла
                    Scalar d = val - (rhs.ri - act1.le);
                    if (StrongGt(d, 0) and (d < val)) {
                        info.UpdCoeff(ic, iv, val - d);
                        n_reduced++;
                    }
                }
            } else {
                // if (rhs.ri == kInf) {
                //     Scalar d = val - (rhs.ri - act1.le);
                //     if (StrongGt(d, 0)) {
                //         act.ri -= d;
                //         el.valueRef() -= d;
                //         model.GetAc().GetCol(iv).AtRef(ic) -= d;
                //         n_reduced++;
                //     }
                // }
            }
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
