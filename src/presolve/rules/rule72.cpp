#include <utility>

#include "reshala/presolve/rules.h"

namespace reshala {

struct Bounder {
    Bounder(const ModelTracker& t)
        : tracker(t),
          activities(t.GetActivities()),
          var_bounds(t.GetModel().GetDomain().GetBounds()) {}

    std::vector<Activity> activities;
    std::vector<Bounds> var_bounds;
    const ModelTracker& tracker;

    bool Propagate(Index iv, const Bounds& new_bnd) {
        const MilpModel& model = tracker.GetModel();
        const Bounds old_bnd = var_bounds[iv];
        var_bounds[iv] = new_bnd;

        std::vector<Index> changed_cons;

        for (SvIterator el(model.GetCol(iv)); el; ++el) {
            Index ic = el.index();
            Scalar val = el.value();
            if (tracker.GetConMask(ic)) continue;

            Activity act = activities[ic];
            Bounds old_range = act.GetRange();
            act.RmTerm(val, old_bnd);
            act.AddTerm(val, new_bnd);
            Bounds range = act.GetRange();
            if (StrongGt(range.le, range.ri)) {
                return false;
            }

            if (range.le != old_range.le or range.ri != old_range.ri) {  // todo weak comparison
                activities[ic] = act;
                changed_cons.push_back(ic);
            }
        }

        for (auto ic : changed_cons) {
            for (SvIterator el(model.GetRow(ic)); el; ++el) {
                Index iv1 = el.index();
                if (tracker.GetVarMask(iv1)) continue;

                Scalar val = el.value();
                const Bounds& bnd = var_bounds[iv1];
                Bounds derived = tracker.DeriveBounds(ic, iv1, activities[ic], bnd, val);

                if (StrongGt(derived.le, bnd.le) or StrongLt(derived.ri, bnd.ri)) {
                    Bounds new_bnd = {std::max(bnd.le, derived.le), std::min(bnd.ri, derived.ri)};

                    if (new_bnd.le > new_bnd.ri) {
                        return false;
                    }
                    var_bounds[iv1] = new_bnd;
                }
            }
        }
        return true;
    }
};

RuleResult Rule72::Apply(ModelTracker& tracker) {
    const MilpModel& model = tracker.GetModel();
    Index n_reduced = 0;

    for (Index iv = 0; iv < model.GetNVars(); iv++) {
        if (tracker.GetVarMask(iv)) continue;
        if (!model.IsBinary(iv)) continue;

        std::array<Bounder, 2> bounders = {Bounder{tracker}, Bounder{tracker}};
        std::array<bool, 2> results;

        for (Index i = 0; i < 2; i++) {
            results[i] = bounders[i].Propagate(iv, {Scalar(i), Scalar(i)});
        }

        if (!results[0] and !results[1]) {
            // printf("Infeasible probing by x%d\n", iv);
            return RuleResult::kInfeasible;
        }
        if (!results[0]) {
            // printf("x%d -> 0\n", iv);
            tracker.UpdVarBounds(iv, {Scalar(0), Scalar(0)});
            n_reduced++;
            continue;
        }
        if (!results[1]) {
            // printf("x%d -> 1\n", iv);
            tracker.UpdVarBounds(iv, {Scalar(1), Scalar(1)});
            n_reduced++;
            continue;
        }

        for (Index iv1 = 0; iv1 < model.GetNVars(); iv1++) {
            const Bounds& bnd = model.GetBounds(iv1);
            Bounds derived = {
                std::min(bounders[0].var_bounds[iv1].le, bounders[1].var_bounds[iv1].le),
                std::max(bounders[0].var_bounds[iv1].ri, bounders[1].var_bounds[iv1].ri),
            };

            if (StrongGt(derived.le, bnd.le) or StrongLt(derived.ri, bnd.ri)) {
                Bounds new_bnd = {std::max(bnd.le, derived.le), std::min(bnd.ri, derived.ri)};
                if (StrongGt(new_bnd.le, new_bnd.ri)) {
                    return RuleResult::kInfeasible;
                }
                tracker.UpdVarBounds(iv1, derived);
                n_reduced++;
            }
        }
    }

    return n_reduced > 0 ? RuleResult::kReduced : RuleResult::kUnchanged;
}

}  // namespace reshala
