#include "reshala/numerics.h"
#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule35::Apply(ModelTracker& tracker) {
    const MilpModel& model = tracker.GetModel();
    Index n_reduced = 0;

    auto obj_gcd = GetGcd(model.GetObj().coefficients);
    if (obj_gcd > 1) {
        tracker.ScaleObj(1. / obj_gcd);
    } else if (!model.ObjIsInteger()) {
        Index min_exp, max_exp;
        Index scale = 0;
        GetExpRange(model.GetObj().coefficients, min_exp, max_exp);
        scale = std::max(scale, min_exp);
        scale = std::min(scale, max_exp);
        if (scale != 0) {
            tracker.ScaleObjExp(-scale);
        }
    }

    for (Index ic = 0; ic < model.GetNCons(); ic++) {
        if (tracker.GetConMask(ic)) continue;

        const SparseVector& row = model.GetRow(ic);
        auto gcd = GetGcd(row.values());
        if (gcd > 1) {
            tracker.ScaleRow(ic, 1. / gcd);
            n_reduced++;
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
