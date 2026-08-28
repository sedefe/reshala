#include "reshala/heuristics/diving.h"

#include <array>

#include "reshala/lp/dual_simplex.h"

namespace reshala {

Solution Diving::InternalRun(const MilpModel& model, const Solution& relaxed,
                             const MipTracker& mip_tracker) {
    if (relaxed.y >= mip_tracker.GetCutoff()) {
        return InfeasibleSolution();
    }

    Solution sol;

    MilpModel model_copy = model;

    auto n_fixed = Fixing(fixing_type_, model_copy, relaxed.x);

    Presolver presolver(model_copy);
    if (n_fixed > 0) {
        LpStatus presolve_status = presolver.Presolve(false, RuleType::kFast);
        if (presolve_status != LpStatus::kUnknown) {
            return presolver.Postsolve({presolve_status, {}, {}});
        }
    }

    DualSimplex ds;
    ds.SetModel(model_copy);
    sol = ds.Solve(false, mip_tracker.GetCutoff());

    while (true) {
        if (sol.status != LpStatus::kOptimal) break;
        if (sol.y >= mip_tracker.GetCutoff()) {
            sol.status = LpStatus::kDropped;
            break;
        }
        if (model_copy.IsIntegerFeasible(sol.x)) break;

        Index cand = GetCandidate(model_copy, relaxed, sol);

        Scalar lb = Floor(sol.x[cand]);
        Scalar rb = lb + 1;

        Bounds bnd = model_copy.GetBounds(cand);
        if (sol.x[cand] - bnd.le < bnd.ri - sol.x[cand])
            bnd.ri = lb;
        else
            bnd.le = rb;
        ds.SetBounds(cand, bnd);
        sol = ds.Solve(true, mip_tracker.GetCutoff());
    }

    return presolver.Postsolve(sol);
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
