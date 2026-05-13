#include "reshala/milp/milp.h"

namespace reshala {

MilpSolver::MilpSolver(MilpModel& model)
    : model(model), presolver(model), mip_state(model), bnb(model, mip_state) {}

Solution MilpSolver::Solve() {
    presolver.Presolve();

    DualSimplex ds(model);
    Solution sol = ds.Solve();
    std::cout << "Root LP: " << sol.y << "\n";


    // Wrap into MipState
    if (sol.status != LpStatus::kOptimal) {
        return sol;
    }

    if (model.IsIntegerFeasible(sol.x)) {
        return presolver.Postsolve(sol);
    }

    bnb.Solve(sol);

    if (mip_state.GetBestSol().status == LpStatus::kOptimal) {
        return presolver.Postsolve(mip_state.GetBestSol());
    } else {
        return {LpStatus::kInfeasible};
    }
}

}  // namespace reshala
