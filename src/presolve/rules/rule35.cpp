#include <limits>
#include <numeric>

#include "reshala/presolve/rules.h"

namespace reshala {

Index GetGcd(const std::vector<Scalar>& vec) {
    Index gcd = 0;
    for (auto x : vec) {
        if ((std::abs(x) <= kMaxInt) and GetFraction(x) == 0.0) {
            gcd = std::gcd(gcd, Index(x));
        } else {
            return 0;
        }
    }
    return gcd;
}

void GetExpRange(const std::vector<Scalar>& vec, Index& min_exp, Index& max_exp) {
    min_exp = std::numeric_limits<Index>::max();
    max_exp = std::numeric_limits<Index>::min();
    Index exp;

    for (auto x : vec) {
        if (x == 0) continue;
        std::frexp(x, &exp);
        min_exp = std::min(min_exp, exp);
        max_exp = std::max(max_exp, exp);
    }
}

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
