#include "reshala/heuristics/diving.h"

#include <array>

#include "reshala/lp/dual_simplex.h"

namespace reshala {

Solution Diving::InternalRun(const Solution& relaxed) {
    if (relaxed.y >= ctx_.GetMip().GetCutoff()) {
        return Solution::Infeasible();
    }

    Solution sol;

    MilpModel model = ctx_.GetModel();

    auto n_fixed = Fixing(fixing_type_, model, relaxed.x);

    Presolver presolver(model);
    bool do_presolve = n_fixed > 0;
    if (do_presolve) {
        LpStatus presolve_status = presolver.Presolve(false, RuleType::kFast);
        if (presolve_status != LpStatus::kUnknown) {
            return presolver.Postsolve({presolve_status, {}, {}});
        }
    }

    DualSimplex ds(ctx_);
    ds.SetModel(model);
    sol = ds.Solve(false);

    while (true) {
        if (sol.status != LpStatus::kOptimal) break;
        if (sol.y >= ctx_.GetMip().GetCutoff()) {
            sol.status = LpStatus::kDropped;
            break;
        }
        if (model.IsIntegerFeasible(sol.x)) break;

        Index cand = GetCandidate(model, relaxed, sol);

        Scalar lb = Floor(sol.x[cand]);
        Scalar rb = lb + 1;

        Bounds bnd = model.GetBounds(cand);
        if (sol.x[cand] - bnd.le < bnd.ri - sol.x[cand])
            bnd.ri = lb;
        else
            bnd.le = rb;
        ds.SetBounds(cand, bnd);
        sol = ds.Solve(true);
    }

    return do_presolve ? presolver.Postsolve(sol) : sol;
}

Index Diving::GetCandidate(const MilpModel& model, const Solution& relaxed, const Solution& sol) {
    // Todo: enhance
    Index candidate = -1;
    Scalar min_fraction = kInf;

    for (Index iv = 0; iv < sol.x.size(); ++iv) {
        if (!model.GetIntegrality(iv)) continue;
        Scalar current_fraction = MinFraction(sol.x[iv]);
        if (IsZero(current_fraction)) continue;

        if (current_fraction < min_fraction) {
            min_fraction = current_fraction;
            candidate = iv;
        }
    }

    return candidate;
}

}  // namespace reshala
