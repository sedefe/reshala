#include "reshala/milp/milp.h"

namespace reshala {

MilpSolver::MilpSolver(Reshala& ctx) : ctx_(ctx) {}

Solution MilpSolver::Solve() {
    MilpModel& model = ctx_.GetModel();
    Presolver& presolver = ctx_.GetPresolver();
    DualSimplex& ds = ctx_.GetDs();
    MipTracker& mip_tracker = ctx_.GetMip();
    mip_tracker.Init();
    HeuristicManager& heur_manager = ctx_.GetHeurMng();

    auto [presolve_status, t_presolve] =
        MEASURE_TIME(presolver.Presolve(true, RuleType::kExhaustive));
    std::cout << "Presolve finished in " << t_presolve << " ms\n";
    if (presolve_status != LpStatus::kUnknown) {
        return presolver.Postsolve({presolve_status, {}, {}});
    }

    ds.SetModel(model);
    auto [sol, t_root] = MEASURE_TIME(ds.Solve(false));
    std::cout << "Root LP: " << sol.y << ", " << t_root << " ms, " << ds.GetStats().n_iter
              << " iterations\n";

    mip_tracker.TestPrimal(sol);
    mip_tracker.UpdDual(sol.y);
    if (mip_tracker.Converged()) {
        return presolver.Postsolve(mip_tracker.GetBestSol());
    }

    heur_manager.Run(HeuristicTrigger::kRoot, sol);
    if (mip_tracker.Converged()) {
        return presolver.Postsolve(mip_tracker.GetBestSol());
    }

    Cutter& cutter = ctx_.GetCutter();
    cutter.Run(sol);
    if (mip_tracker.Converged()) {
        return presolver.Postsolve(mip_tracker.GetBestSol());
    }

    BnbSolver& bnb = ctx_.GetBnb();
    bnb.Solve(sol);

    return presolver.Postsolve(mip_tracker.GetBestSol());
}

}  // namespace reshala
