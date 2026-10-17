#include "reshala/numerics.h"
#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule35::Apply(ModelTracker& tracker, Index ic) {
    const Index kGcdDenominator = 600;
    const Index kFpExpScalingLim = 5;

    const MilpModel& model = tracker.GetModel();

    // Todo run once
    {  // objective
        auto obj_gcd = GetGcd(model.GetObj().coefficients, kGcdDenominator);
        if (obj_gcd != 0 and obj_gcd != kGcdDenominator) {
            tracker.ScaleObj(Scalar(kGcdDenominator) / obj_gcd);  // Todo rational scaling
        }
    }

    {  // constraints
        const SparseVector& row = model.GetRow(ic);
        auto gcd = GetGcd(row.values(), kGcdDenominator);
        if (gcd != 0 and gcd != kGcdDenominator) {
            tracker.ScaleCon(ic, Scalar(kGcdDenominator) / gcd);
            return RuleResult::kReduced;
        }
    }

    return RuleResult::kUnchanged;
}

}  // namespace reshala
