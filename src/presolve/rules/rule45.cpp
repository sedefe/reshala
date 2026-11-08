#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule45::Apply(ModelTracker& tracker) {
    // Todo: avoid CalcActivities()

    const MilpModel& model = tracker.GetModel();
    Index n_reduced = 0;

    for (Index iv = 0; iv < model.GetNVars(); iv++) {
        if (tracker.VarIsDeleted(iv)) continue;
        bool var_integral = model.GetIntegrality(iv);

        Bounds implied = {-kInf, kInf};
        const auto& bnd = model.GetBounds(iv);
        const auto& col = model.GetCol(iv);

        for (SvIterator el(col); el; ++el) {
            implied = BoundsIntersection(
                implied, tracker.DeriveBounds(el.index(), iv, tracker.GetActivity(el.index()), bnd,
                                              el.value()));
        }

        if (WeakGe(implied.le, bnd.le) and WeakLe(implied.ri, bnd.ri)) {
            for (SvIterator el(col); el; ++el) {
                Index ic = el.index();
                if (tracker.ConIsDeleted(ic)) continue;

                const auto& row = model.GetRow(ic);
                const auto& rhs = model.GetRhs(ic);
                if (!WeakEq(rhs.le, rhs.ri)) continue;

                Scalar val = el.value();
                if (std::abs(val) < kMinAbsPivot) continue;

                bool eligible = true;
                SparseVector sv(model.GetNVars());
                for (SvIterator el(row); el and eligible; ++el) {
                    if (tracker.VarIsDeleted(el.index())) continue;

                    Scalar coef = -el.value() / val;
                    if (el.index() != iv) sv.Push(el.index(), coef);

                    if (var_integral) {
                        eligible &= IsZero(MinFraction(coef));
                    }
                }
                Scalar a = (rhs.le + rhs.ri) / 2 / val;

                {  // integrality check
                    eligible &= IsZero(MinFraction(a));
                    if (!eligible) continue;
                }

                Index fillin = sv.Size() * (col.Size() - 1);
                if (fillin > model.GetNnz() * kMaxFillin) continue;

                tracker.MaskCon(ic);
                tracker.LinCombSub(iv, sv, a);
                tracker.MaskVar(iv);
                tracker.RebuildAc();
                tracker.CalcActivities();
                for (SvIterator el(row); el; ++el) {
                    if (tracker.VarIsDeleted(el.index())) continue;
                    tracker.UpdLocks(el.index());
                }

                n_reduced++;
                break;
            }
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
