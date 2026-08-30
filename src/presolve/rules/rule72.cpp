#include <array>

#include "reshala/presolve/rules.h"

namespace reshala {

struct Bounder {
    Bounder(const ModelTracker& t)
        : tracker(t), con_mask(t.GetModel().GetNCons()), var_mask(t.GetModel().GetNVars()) {
        changed_cons.reserve(t.GetModel().GetNCons());
        changed_vars.reserve(t.GetModel().GetNVars());
    }
    void Reset() {
        activities = tracker.GetActivities();
        var_bounds = tracker.GetModel().GetDomain().GetBounds();
    }

    const ModelTracker& tracker;

    std::vector<Activity> activities;
    std::vector<Bounds> var_bounds;

    std::vector<Index> changed_cons;
    BitMask con_mask;
    std::vector<Index> changed_vars;
    BitMask var_mask;

    bool Propagate(Index iv_start, const Bounds& new_bnd) {
        const MilpModel& model = tracker.GetModel();

        const Index kMaxIters = 5;
        Index n_iter = 0;

        changed_vars.clear();
        var_mask.Clear();
        changed_vars.push_back(iv_start);

        std::vector<Bounds> old_bounds = var_bounds;
        var_bounds[iv_start] = new_bnd;

        while (!changed_vars.empty() and n_iter < kMaxIters) {
            n_iter++;

            changed_cons.clear();
            con_mask.Clear();
            for (auto iv : changed_vars) {  // var_bounds -> activities
                for (SvIterator el(model.GetCol(iv)); el; ++el) {
                    Index ic = el.index();
                    Scalar val = el.value();
                    if (tracker.GetConMask(ic)) continue;

                    Activity act = activities[ic];
                    Bounds old_range = act.GetRange();
                    act.RmTerm(val, old_bounds[iv]);
                    act.AddTerm(val, var_bounds[iv]);
                    Bounds new_range = act.GetRange();

                    if (StrongGt(new_range.le, old_range.le) or
                        StrongLt(new_range.ri, old_range.ri)) {
                        activities[ic] = act;

                        if (StrongGt(new_range.le, model.GetRhs(ic).ri) or
                            StrongLt(new_range.ri, model.GetRhs(ic).le)) {
                            return false;
                        }

                        if (!con_mask.Get(ic)) {
                            con_mask.Set(ic);
                            changed_cons.push_back(ic);
                        }
                    }
                }
            }

            old_bounds = var_bounds;

            changed_vars.clear();
            var_mask.Clear();
            for (auto ic : changed_cons) {  // activities -> var_bounds
                for (SvIterator el(model.GetRow(ic)); el; ++el) {
                    Index iv = el.index();
                    if (tracker.GetVarMask(iv)) continue;

                    Scalar val = el.value();
                    const Bounds& old_bnd = old_bounds[iv];
                    Bounds& curr_bnd = var_bounds[iv];
                    Bounds derived = tracker.DeriveBounds(ic, iv, activities[ic], old_bnd, val);

                    if (StrongGt(derived.le, curr_bnd.le) or StrongLt(derived.ri, curr_bnd.ri)) {
                        curr_bnd = {std::max(curr_bnd.le, derived.le),
                                    std::min(curr_bnd.ri, derived.ri)};
                        if (StrongGt(curr_bnd.le, curr_bnd.ri)) {
                            return false;
                        }
                        if (!var_mask.Get(iv)) {
                            var_mask.Set(iv);
                            changed_vars.push_back(iv);
                        }
                    }
                }
            }
        }
        return true;
    }
};

RuleResult Rule72::Apply(ModelTracker& tracker) {
    const MilpModel& model = tracker.GetModel();
    Index n_reduced = 0;
    tracker.GetImplications().clear();

    std::array<Bounder, 2> bounders = {Bounder{tracker}, Bounder{tracker}};
    std::array<bool, 2> results;

    for (Index iv = 0; iv < model.GetNVars(); iv++) {
        if (tracker.GetVarMask(iv)) continue;
        if (!model.IsBinary(iv)) continue;

        for (Index i = 0; i < 2; i++) {
            bounders[i].Reset();
            results[i] = bounders[i].Propagate(iv, {Scalar(i), Scalar(i)});
        }

        if (!results[0] and !results[1]) {
            return RuleResult::kInfeasible;
        }
        if (!results[0]) {
            tracker.UpdVarBounds(iv, {Scalar(1), Scalar(1)});
            n_reduced++;
            continue;
        }
        if (!results[1]) {
            tracker.UpdVarBounds(iv, {Scalar(0), Scalar(0)});
            n_reduced++;
            continue;
        }

        for (Index iv1 = 0; iv1 < model.GetNVars(); iv1++) {
            if (tracker.GetVarMask(iv1)) continue;
            if (iv == iv1) continue;
            const Bounds& bnd = model.GetBounds(iv1);
            const Bounds& bnd0 = bounders[0].var_bounds[iv1];
            const Bounds& bnd1 = bounders[1].var_bounds[iv1];
            Bounds derived = {
                std::min(bnd0.le, bnd1.le),
                std::max(bnd0.ri, bnd1.ri),
            };

            // Check if we can fix or substitute x <- a*y + b
            if (WeakEq(bnd0.le, bnd0.ri) and WeakEq(bnd1.le, bnd1.ri)) {
                Scalar y0 = (bnd0.le + bnd0.ri) / 2;
                Scalar y1 = (bnd1.le + bnd1.ri) / 2;
                if (WeakEq(y0, y1)) {
                    tracker.FixVar(iv1, (y0 + y1) / 2);
                } else {
                    tracker.SimpleSub(iv1, y1 - y0, iv, y0);
                }
                n_reduced++;
                continue;
            }

            // Check if we can strengthen the bounds
            if (StrongGt(derived.le, bnd.le) or StrongLt(derived.ri, bnd.ri)) {
                Bounds new_bnd = {std::max(bnd.le, derived.le), std::min(bnd.ri, derived.ri)};
                if (StrongGt(new_bnd.le, new_bnd.ri)) {
                    return RuleResult::kInfeasible;
                }
                tracker.UpdVarBounds(iv1, std::move(new_bnd));
                n_reduced++;
                continue;
            }

            if (StrongGt(bnd0.le, bnd.le)) {  // x = 0 => y >= b
                tracker.GetImplications().push_back({iv, false, iv1, true, bnd0.le});
            }
            if (StrongLt(bnd0.ri, bnd.ri)) {  // x = 0 => y <= b
                tracker.GetImplications().push_back({iv, false, iv1, false, bnd0.ri});
            }
            if (StrongGt(bnd1.le, bnd.le)) {  // x = 1 => y >= b
                tracker.GetImplications().push_back({iv, true, iv1, true, bnd1.le});
            }
            if (StrongLt(bnd1.ri, bnd.ri)) {  // x = 1 => y <= b
                tracker.GetImplications().push_back({iv, true, iv1, false, bnd1.ri});
            }
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
