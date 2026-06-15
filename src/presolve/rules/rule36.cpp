#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule36::Apply(ModelTracker& tracker) {
    const MilpModel& model = tracker.GetModel();
    Index n_reduced = 0;
    for (Index ic = 0; ic < model.GetNCons(); ic++) {
        if (tracker.GetConMask(ic)) continue;

        auto rhs = model.GetRhs(ic);
        if (!IsZero(rhs.ri - rhs.le)) continue;
        Scalar b = (rhs.ri + rhs.le) / 2;
        auto range = tracker.GetActivity(ic).GetRange();
        if (IsZero(range.le + range.ri - 2 * b)) {
            const auto& row = model.GetRow(ic);
            Index iv_base = -1;
            Scalar a;
            for (SvIterator el(row); el; ++el) {
                if (tracker.GetVarMask(el.index())) continue;
                if (model.IsBinary(el.index())) {
                    if (IsZero(std::abs(el.value()) - (range.ri - b))) {
                        iv_base = el.index();
                        a = el.value();
                        break;
                    }
                }
            }
            if (iv_base >= 0) {
                tracker.MaskCon(ic);  // Сразу маскируем, чтобы при подстановке не менять эту строку
                for (SvIterator el(row); el; ++el) {
                    if (tracker.GetVarMask(el.index())) continue;
                    if (el.index() == iv_base) continue;

                    const auto& bnd = model.GetBounds(el.index());
                    Scalar range = bnd.ri - bnd.le;
                    if (a >= 0) {
                        if (el.value() >= 0) {  // x + y = 1 => y = -x + 1
                            if (!tracker.SimpleSub(el.index(), -range, iv_base, bnd.ri))
                                return RuleResult::kInfeasible;
                        } else {  // x - y = 0 => y = x + 0
                            if (!tracker.SimpleSub(el.index(), range, iv_base, bnd.le))
                                return RuleResult::kInfeasible;
                        }
                    } else {
                        if (el.value() >= 0) {  // -x + y = 0 => y = x
                            if (!tracker.SimpleSub(el.index(), range, iv_base, bnd.le))
                                return RuleResult::kInfeasible;
                        } else {  // -x - y = -1 => y = 1-x
                            if (!tracker.SimpleSub(el.index(), -range, iv_base, bnd.ri))
                                return RuleResult::kInfeasible;
                        }
                    }
                }
                n_reduced++;
                continue;
            }
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
