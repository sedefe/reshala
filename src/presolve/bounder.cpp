#include "reshala/presolve/bounder.h"

#include "reshala/presolve/model_tracker.h"

namespace reshala {

Bounder::Bounder(const ModelTracker& t)
    : tracker(t),
      changed_con_mask(t.GetModel().GetNCons()),
      changed_var_mask(t.GetModel().GetNVars()),
      all_changed_con_mask(t.GetModel().GetNCons()),
      all_changed_var_mask(t.GetModel().GetNVars()),
      redundant_con_mask(t.GetModel().GetNCons()) {
    changed_cons.reserve(t.GetModel().GetNCons());
    changed_vars.reserve(t.GetModel().GetNVars());
}

void Bounder::Reset() {
    activities = tracker.GetActivities();
    domain = tracker.GetModel().GetDomain();
    redundant_con_mask = tracker.GetConMask();

    changed_vars.clear();
    changed_var_mask = tracker.GetVarMask();
    changed_cons.clear();
    changed_con_mask = tracker.GetConMask();

    all_changed_vars.clear();
    all_changed_var_mask = tracker.GetVarMask();
    all_changed_cons.clear();
    all_changed_con_mask = tracker.GetConMask();
}

bool Bounder::Propagate(Index iv_start, const Bounds& new_bnd) {
    const MilpModel& model = tracker.GetModel();

    const Index kMaxIters = 15;
    Index n_iter = 0;

    changed_vars.push_back(iv_start);

    Domain old_domain = domain;
    domain.SetBounds(iv_start, new_bnd);

    while (!changed_vars.empty() and n_iter < kMaxIters) {
        n_iter++;

        changed_cons.clear();
        changed_con_mask = tracker.GetConMask();
        for (auto iv : changed_vars) {  // domain -> activities
            for (SvIterator el(model.GetCol(iv)); el; ++el) {
                Index ic = el.index();
                Scalar val = el.value();
                if (redundant_con_mask.Get(ic)) continue;

                Activity act = activities[ic];
                Bounds old_range = act.GetRange();
                act.RmTerm(val, old_domain.GetBounds(iv));
                act.AddTerm(val, domain.GetBounds(iv));
                Bounds new_range = act.GetRange();

                if (WeakGe(new_range.le, model.GetRhs(ic).le) and
                    WeakLe(new_range.ri, model.GetRhs(ic).ri)) {
                    redundant_con_mask.Set(ic);
                    if (!all_changed_con_mask.Get(ic)) {
                        all_changed_con_mask.Set(ic);
                        all_changed_cons.push_back(ic);
                    }
                    continue;
                }

                if (StrongGt(new_range.le, old_range.le) or StrongLt(new_range.ri, old_range.ri)) {
                    activities[ic] = act;

                    if (StrongGt(new_range.le, model.GetRhs(ic).ri) or
                        StrongLt(new_range.ri, model.GetRhs(ic).le)) {
                        return false;
                    }

                    if (!changed_con_mask.Get(ic)) {
                        changed_con_mask.Set(ic);
                        changed_cons.push_back(ic);
                        if (!all_changed_con_mask.Get(ic)) {
                            all_changed_con_mask.Set(ic);
                            all_changed_cons.push_back(ic);
                        }
                    }
                }
            }
        }

        old_domain = domain;

        changed_vars.clear();
        changed_var_mask.Clear();
        for (auto ic : changed_cons) {  // activities -> domain
            for (SvIterator el(model.GetRow(ic)); el; ++el) {
                Index iv = el.index();
                if (tracker.GetVarMask(iv)) continue;

                Scalar val = el.value();
                const Bounds& old_bnd = old_domain.GetBounds(iv);
                Bounds curr_bnd = domain.GetBounds(iv);
                Bounds derived = tracker.DeriveBounds(ic, iv, activities[ic], old_bnd, val);

                if (StrongGt(derived.le, curr_bnd.le) or StrongLt(derived.ri, curr_bnd.ri)) {
                    domain.SetBounds(
                        iv, {std::max(curr_bnd.le, derived.le), std::min(curr_bnd.ri, derived.ri)});
                    if (StrongGt(curr_bnd.le, curr_bnd.ri)) {
                        return false;
                    }
                    if (!changed_var_mask.Get(iv)) {
                        changed_var_mask.Set(iv);
                        changed_vars.push_back(iv);
                        if (!all_changed_var_mask.Get(iv)) {
                            all_changed_var_mask.Set(iv);
                            all_changed_vars.push_back(iv);
                        }
                    }
                }
            }
        }
    }
    return true;
}

}  // namespace reshala
