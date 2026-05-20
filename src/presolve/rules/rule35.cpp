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

RuleResult Rule35::Apply(ModelInfo& info, std::vector<std::unique_ptr<Transform>>& transforms) {
    MilpModel& model = info.GetModel();
    Index n_reduced = 0;

    auto obj_gcd = GetGcd(info.GetModel().GetObj().coefficients);
    if (obj_gcd > 1) {
        for (auto& x : info.GetModel().GetObj().coefficients) x /= obj_gcd;
        info.GetModel().GetObj().mult *= obj_gcd;
    }

    for (Index ic = 0; ic < model.GetNCons(); ic++) {
        if (info.GetConMask(ic)) continue;

        SparseVector& row = info.GetModel().GetAr().GetRow(ic);
        auto gcd = GetGcd(row.values());
        if (gcd > 1) {
            row *= (1. / gcd);
            info.GetActivity(ic).le /= gcd;
            info.GetActivity(ic).ri /= gcd;
            info.GetModel().GetRhs(ic).le /= gcd;
            info.GetModel().GetRhs(ic).ri /= gcd;
            n_reduced++;
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
