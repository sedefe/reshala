#include "reshala/heuristics/diving.h"

#include <array>

#include "reshala/lp/dual_simplex.h"

namespace reshala {

Solution Diving::InternalRun(const MilpModel& model, const Solution& relaxed,
                             const MipState& mip_state) {
    Solution sol;

    MilpModel model_copy = model;

    Fixing(type_, model_copy, relaxed.x);

    Presolver presolver(model_copy);
    LpStatus presolve_status = presolver.Presolve(false);
    if (presolve_status != LpStatus::kUnknown) {
        return presolver.Postsolve({presolve_status, {}, {}});
    }

    DualSimplex ds;
    ds.SetModel(model_copy);
    sol = ds.Solve(false);

    while (true) {
        if (sol.status == LpStatus::kInfeasible) break;
        if (model_copy.IsIntegerFeasible(sol.x)) break;

        Index cand = GetCandidate(model_copy, relaxed, sol);

        std::array<Bounds, 2> bounds_priority;
        Scalar lb = Floor(sol.x[cand]);
        Scalar rb = lb + 1;

        // Todo: enhance logic, don't set constant bounds
        const Bounds& bnd = model_copy.GetBounds(cand);
        if (sol.x[cand] - bnd.le > bnd.ri - sol.x[cand])
            bounds_priority = {{{lb, lb}, {rb, rb}}};
        else
            bounds_priority = {{{rb, rb}, {lb, lb}}};

        // Todo: choose the best child
        for (Index i = 0; i < 2; ++i) {
            ds.SetBounds(cand, bounds_priority[i]);
            sol = ds.Solve(true);
            if (sol.status == LpStatus::kOptimal) break;
        }
    }

    return presolver.Postsolve(sol);
}

Index Diving::GetCandidate(const MilpModel& model, const Solution& relaxed, const Solution& sol) {
    // Todo: enhance
    Index candidate = -1;
    Scalar max_fraction = -kInf;

    for (Index iv = 0; iv < sol.x.size(); ++iv) {
        if (!model.GetIntegrality(iv)) continue;
        Scalar current_fraction = GetFraction(sol.x[iv]);
        if (current_fraction > max_fraction) {
            max_fraction = current_fraction;
            candidate = iv;
        }
    }

    return candidate;
}

}  // namespace reshala
