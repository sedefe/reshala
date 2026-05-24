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

RuleResult Rule35::Apply(ModelTracker& info, std::vector<std::unique_ptr<Transform>>& transforms) {
    const MilpModel& model = info.GetModel();
    Index n_reduced = 0;

    auto obj_gcd = GetGcd(model.GetObj().coefficients);
    if (obj_gcd > 1) {
        info.ScaleObj(1. / obj_gcd);
    }

    for (Index ic = 0; ic < model.GetNCons(); ic++) {
        if (info.GetConMask(ic)) continue;

        const SparseVector& row = model.GetAr().GetRow(ic);
        auto gcd = GetGcd(row.values());
        if (gcd > 1) {
            info.ScaleRow(ic, 1. / gcd);
            n_reduced++;
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
