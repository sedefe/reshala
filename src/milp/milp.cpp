#include "reshala/milp/milp.h"

namespace reshala {

MilpSolver::MilpSolver(MilpModel& model)
    : model(model),
      mip_state(model),
      presolver(model),
      heuristics(model, mip_state),
      bnb(model, mip_state) {}

Solution MilpSolver::Solve() {
    LpStatus presolve_status = presolver.Presolve();
    if (presolve_status != LpStatus::kUnknown) {
        return presolver.Postsolve({presolve_status, model.GetObj().c0, {}});
    }

    DualSimplex ds(model);
    auto [sol, duration] = MEASURE_TIME(ds.Solve());
    std::cout << "Root LP: " << sol.y << ", " << duration.count() / 1e3 << " ms\n";

    mip_state.TestPrimal(sol);
    mip_state.UpdDual(sol.y);
    if (mip_state.Converged()) {
        return presolver.Postsolve(mip_state.GetBestSol());
    }

    Solution sol_h = heuristics.Run(sol);
    mip_state.TestPrimal(sol_h);
    if (mip_state.Converged()) {
        return presolver.Postsolve(mip_state.GetBestSol());
    }

    bnb.Solve(sol);

    return presolver.Postsolve(mip_state.GetBestSol());
}

}  // namespace reshala
