#include "reshala/heuristics/diving.h"

#include "reshala/lp/dual_simplex.h"

namespace reshala {

Solution Diving::InternalRun(MilpModel& model, const Solution& relaxed, const MipState& mip_state) {
    Solution res;

    Fixing(type_, model, relaxed.x);

    Presolver presolver(model);
    LpStatus presolve_status = presolver.Presolve();
    if (presolve_status != LpStatus::kUnknown) {
        return presolver.Postsolve({presolve_status, model.GetObj().c0, {}});
    }

    MipState mip_state_copy(model, mip_state.GetBestSol(), mip_state.GetDual());
    DualSimplex ds(model);
    while (true) {
        Solution sol = ds.Solve(false);

        if (sol.status == LpStatus::kInfeasible) break;
        if (model.IsIntegerFeasible(sol.x)) break;
        mip_state_copy.UpdDual(sol.y);
        if (mip_state_copy.Converged()) {
            break;
        }

        Index cand = GetCandidate(model, relaxed, sol);

        std::array<Bounds, 2> bounds_priority;
        Scalar lb = Floor(sol.x[cand]);
        Scalar rb = lb + 1;

        // Todo: enhance
        const Bounds& bnd = model.GetBounds(cand);
        if (sol.x[cand] - bnd.le > bnd.ri - sol.x[cand])
            bounds_priority = {{{lb, lb}, {rb, rb}}};
        else
            bounds_priority = {{{rb, rb}, {lb, lb}}};

        Solution child;
        // Todo: choose the best child
        for (Index i = 0; i < 2; ++i) {
            DualSimplex ds1(model);
            model.SetBounds(cand, bounds_priority[i]);
            child = ds1.Solve(false);
            if (child.status == LpStatus::kOptimal) break;
        }
    }

    return presolver.Postsolve(res);
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
