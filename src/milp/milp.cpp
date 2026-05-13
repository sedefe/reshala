#include "reshala/milp/milp.h"

namespace reshala {

MilpSolver::MilpSolver(MilpModel& model)
    : model(model), mip_state(model), presolver(model), bnb(model, mip_state) {}

Solution MilpSolver::Solve() {
    presolver.Presolve();

    DualSimplex ds(model);
    Solution sol = ds.Solve();
    std::cout << "Root LP: " << sol.y << "\n";

    mip_state.TestPrimal(sol);
    mip_state.UpdDual(sol.y);

    if (!mip_state.Converged()) {
        bnb.Solve(sol);
    }

    return presolver.Postsolve(mip_state.GetBestSol());
}

}  // namespace reshala
