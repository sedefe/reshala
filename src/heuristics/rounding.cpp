#include "reshala/heuristics/rounding.h"

namespace reshala {

Solution Rounding::InternalRun(const MilpModel& model, const Solution& relaxed,
                               const MipTracker& mip_tracker) {
    Solution sol;
    sol.status = LpStatus::kInfeasible;
    sol.x = relaxed.x;

    bool eligible = true;
    for (Index iv = 0; iv < model.GetNVars(); iv++) {
        if (!model.GetIntegrality(iv)) continue;
        if (IsZero(MinFraction(relaxed.x[iv]))) continue;

        if (model.GetNLocks(iv, LockType::kUp) == 0) {
            if (model.GetNLocks(iv, LockType::kDown) == 0) {
                // Todo it's either free or does not appear in constraints
                std::cerr << "Handle me\n";
                exit(0);
            } else {
                sol.x[iv] = Ceil(sol.x[iv]);
            }
        } else {
            if (model.GetNLocks(iv, LockType::kDown) == 0) {
                sol.x[iv] = Floor(sol.x[iv]);
            } else {
                eligible = false;
                break;
            }
        }
    }

    if (eligible) {
        sol.status = LpStatus::kOptimal;
        sol.y = model.GetObj().evaluate(sol.x);
    }

    return sol;
}

}  // namespace reshala
