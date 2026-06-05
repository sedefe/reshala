#include "reshala/milp/milp.h"

namespace reshala {

MilpSolver::MilpSolver(MilpModel& model)
    : model(model),
      mip_state(model),
      presolver(model),
      heuristics(model, mip_state),
      bnb(model, *ds, mip_state) {}

Solution MilpSolver::Solve() {
    auto [presolve_status, t_presolve] = MEASURE_TIME(presolver.Presolve());
    std::cout << "Presolve finished in " << t_presolve << " ms\n";
    if (presolve_status != LpStatus::kUnknown) {
        return presolver.Postsolve({presolve_status, model.GetObj().c0, {}});
    }

    ds.emplace(model);
    auto [sol, t_root] = MEASURE_TIME(ds->Solve(false));
    std::cout << "Root LP: " << sol.y << ", " << t_root << " ms, " << ds->GetStats().n_iter
              << " iterations\n";

    mip_state.TestPrimal(sol);
    mip_state.UpdDual(sol.y);
    if (mip_state.Converged()) {
        return presolver.Postsolve(mip_state.GetBestSol());
    }

    heuristics.Run(sol);
    if (mip_state.Converged()) {
        return presolver.Postsolve(mip_state.GetBestSol());
    }

    bnb.Solve(sol);

    return presolver.Postsolve(mip_state.GetBestSol());
}

}  // namespace reshala
