#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule46::Apply(ModelTracker& tracker, Index ic) {
    const MilpModel& model = tracker.GetModel();

    auto rhs = model.GetRhs(ic);
    if (!WeakEq(rhs.ri, rhs.le)) return RuleResult::kUnchanged;
    Scalar b = (rhs.ri + rhs.le) / 2;

    const auto& row = model.GetRow(ic);
    SparseVector compressed_row(model.GetNVars());

    for (SvIterator el(row); el; ++el) {  // Check if there's exactly 2 vars
        if (tracker.VarIsDeleted(el.index())) continue;
        compressed_row.Push(el.index(), el.value());
        if (compressed_row.Size() > 2) break;
    }
    if (compressed_row.Size() != 2) return RuleResult::kUnchanged;

    Index i1 = compressed_row.indices()[0];
    Index i2 = compressed_row.indices()[1];
    Scalar v1 = compressed_row.values()[0];
    Scalar v2 = compressed_row.values()[1];

    if (std::abs(v1) < std::abs(v2)) {
        std::swap(i1, i2);
        std::swap(v1, v2);
    }
    // Now we have v1 x1 + v2 x2 = b, v1 >= v2
    if (!model.GetIntegrality(i1)) {
        if (!tracker.SimpleSub(i1, -v2 / v1, i2, b / v1)) return RuleResult::kInfeasible;
        tracker.MaskCon(ic);
        return RuleResult::kReduced;
    } else {
        if (!model.GetIntegrality(i2)) {
            if (!IsZero(v2 / v1)) {
                if (!tracker.SimpleSub(i2, -v1 / v2, i1, b / v2)) return RuleResult::kInfeasible;
                tracker.MaskCon(ic);
                return RuleResult::kReduced;
            }
        } else {  // v1 x + x2 = b
            v1 /= v2;
            b /= v2;
            if (IsZero(MinFraction(v1))) {
                if (!IsZero(MinFraction(b))) return RuleResult::kInfeasible;
                if (!tracker.SimpleSub(i2, -v1, i1, b)) return RuleResult::kInfeasible;
                tracker.MaskCon(ic);
                return RuleResult::kReduced;
            }
        }
    }

    return RuleResult::kUnchanged;
}

}  // namespace reshala
