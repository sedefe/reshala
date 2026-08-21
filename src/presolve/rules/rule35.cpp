#include "reshala/numerics.h"
#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule35::Apply(ModelTracker& tracker) {
    const Index kGcdDenominator = 600;
    const Index kFpExpScalingLim = 5;

    const MilpModel& model = tracker.GetModel();
    Index n_reduced = 0;

    {  // objective
        auto obj_gcd = GetGcd(model.GetObj().coefficients, kGcdDenominator);
        if (obj_gcd != 0 and obj_gcd != kGcdDenominator) {
            tracker.ScaleObj(Scalar(kGcdDenominator) / obj_gcd);  // Todo rational scaling
        }
    }

    {  // constraints
        for (Index ic = 0; ic < model.GetNCons(); ic++) {
            if (tracker.GetConMask(ic)) continue;

            const SparseVector& row = model.GetRow(ic);
            auto gcd = GetGcd(row.values());
            if (gcd > 1) {
                tracker.ScaleRow(ic, 1. / gcd);
                n_reduced++;
            }
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
