#include <numeric>

#include "reshala/presolve/rules.h"

namespace reshala {

Index GetGcd(const std::vector<Scalar>& vec) {
    Index gcd = 0;
    for (auto x : vec) {
        if (GetFraction(x) == 0.0) {
            gcd = std::gcd(gcd, Index(x));
        } else {
            return 0;
        }
    }
    return gcd;
}

RuleResult Rule35::Apply(ModelTracker& tracker) {
    const MilpModel& model = tracker.GetModel();
    Index n_reduced = 0;

    auto obj_gcd = GetGcd(model.GetObj().coefficients);
    if (obj_gcd > 1) {
        tracker.ScaleObj(1. / obj_gcd);
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
