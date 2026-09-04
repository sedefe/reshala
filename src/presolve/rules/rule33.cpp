#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule33::Apply(ModelTracker& tracker) {
    const MilpModel& model = tracker.GetModel();
    Index n_reduced = 0;

    for (Index ic = 0; ic < model.GetNCons(); ic++) {
        if (tracker.GetConMask(ic)) continue;
        const Bounds& rhs = model.GetRhs(ic);
        if (std::isinf(rhs.le) == std::isinf(rhs.ri)) continue;
        const Activity& act = tracker.GetActivity(ic);

        for (SvIterator el(model.GetRow(ic)); el; ++el) {
            Index iv = el.index();
            if (tracker.GetVarMask(iv)) continue;
            if (!model.IsBinary(iv)) continue;

            const Bounds& bnd = model.GetBounds(iv);
            Scalar val = el.value();
            Activity act1 = act;
            act1.RmTerm(val, bnd);
            auto lhs = act1.GetRange();

            Scalar d = 0;
            if (rhs.le == -kInf) {  // Todo вынести из цикла
                if (val > 0) {
                    d = (rhs.ri - lhs.ri);
                    Scalar new_val = val - d;
                    if (val >= d and StrongGt(d, 0) and StrongGt(val / new_val, 1.)) {
                        tracker.UpdCoeff(ic, iv, new_val);
                        tracker.UpdRhs(ic, {-kInf, lhs.ri});
                        n_reduced++;
                    }
                } else {
                    d = (rhs.ri - val - lhs.ri);
                    Scalar new_val = val + d;
                    if (-val >= d and StrongGt(d, 0) and StrongGt(-val / -new_val, 1.)) {
                        tracker.UpdCoeff(ic, iv, new_val);
                        n_reduced++;
                    }
                }
            } else {
                if (val > 0) {
                    d = (lhs.le + val - rhs.le);
                    Scalar new_val = val - d;
                    if (val >= d and StrongGt(d, 0) and StrongGt(val / new_val, 1.)) {
                        tracker.UpdCoeff(ic, iv, new_val);
                        n_reduced++;
                    }
                } else {
                    d = (lhs.le - rhs.le);
                    Scalar new_val = val + d;
                    if (-val >= d and StrongGt(d, 0) and StrongGt(-val / -new_val, 1.)) {
                        tracker.UpdCoeff(ic, iv, new_val);
                        tracker.UpdRhs(ic, {lhs.le, kInf});
                        n_reduced++;
                    }
                }
            }
        }
    }
    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
