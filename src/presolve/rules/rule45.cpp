#include "reshala/presolve/rules.h"

namespace reshala {

RuleResult Rule45::Apply(ModelTracker& tracker) {
    // Todo: check fillin, avoid CalcActivities()

    const MilpModel& model = tracker.GetModel();
    Index n_reduced = 0;
    Index size = 0;

    for (Index iv = 0; iv < model.GetNVars(); iv++) {
        if (tracker.VarIsDeleted(iv)) continue;
        if (model.GetIntegrality(iv)) continue;

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

                {
                    Index val = el.value();

                    SparseVector sv(model.GetNVars());
                    for (SvIterator el(row); el; ++el) {
                        if (tracker.VarIsDeleted(el.index())) continue;
                        if (el.index() != iv) sv.Push(el.index(), -el.value() / val);
                    }
                    Scalar a = (rhs.le + rhs.ri) / 2 / val;

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
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
