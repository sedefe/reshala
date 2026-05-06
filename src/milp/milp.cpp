#include "reshala/milp/milp.h"

namespace reshala {

MilpSolver::MilpSolver(MilpModel& model)
    : model(model),
      presolver(model),
      ds(model),
      mip_state(-kInf, kInf, model),
      bnb(model, mip_state) {}

Solution MilpSolver::Solve() {
    presolver.Presolve();

    Solution sol = ds.Solve();

    if (sol.status != LpStatus::kOptimal) {
        return sol;
    }

    if (model.IsIntegerFeasible(sol.x)) {
        return sol;
    }

    return bnb.Solve(sol);
}

}  // namespace reshala
