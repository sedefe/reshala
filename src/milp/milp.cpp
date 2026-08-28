#include "reshala/milp/milp.h"

namespace reshala {

MilpSolver::MilpSolver(MilpModel& model)
    : model(model),
      mip_tracker(model),
      presolver(model),
      heur_manager(mip_tracker),
      cutter(model, presolver, ds, mip_tracker, heur_manager),
      bnb(model, ds, mip_tracker, heur_manager) {}

Solution MilpSolver::Solve() {
    auto [presolve_status, t_presolve] =
        MEASURE_TIME(presolver.Presolve(true, RuleType::kExhaustive));
    std::cout << "Presolve finished in " << t_presolve << " ms\n";
    if (presolve_status != LpStatus::kUnknown) {
        return presolver.Postsolve({presolve_status, {}, {}});
    }

    ds.SetModel(model);
    auto [sol, t_root] = MEASURE_TIME(ds.Solve(false, kInf));
    std::cout << "Root LP: " << sol.y << ", " << t_root << " ms, " << ds.GetStats().n_iter
              << " iterations\n";

    mip_tracker.TestPrimal(sol);
    mip_tracker.UpdDual(sol.y);
    if (mip_tracker.Converged()) {
        return presolver.Postsolve(mip_tracker.GetBestSol());
    }

    heur_manager.Run(HeuristicTrigger::kRoot, model, sol);
    if (mip_tracker.Converged()) {
        return presolver.Postsolve(mip_tracker.GetBestSol());
    }

    cutter.Run(sol);
    if (mip_tracker.Converged()) {
        return presolver.Postsolve(mip_tracker.GetBestSol());
    }

    bnb.Solve(sol);

    return presolver.Postsolve(mip_tracker.GetBestSol());
}

void MilpSolver::PrintStats(std::ostream& os) const {
    os << "=== Stats ===\n";
    os << ds.GetLina().GetStats();
    os << ds.GetStats();
    os << ds.GetScaling().stats;
    heur_manager.PrintStats(os);
    os << cutter.GetStats();
    os << bnb.GetStats();
}

}  // namespace reshala
