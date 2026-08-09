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

        Index n_up_locks = 0, n_down_locks = 0;  // Todo keep locks in the model
        for (SvIterator el(model.GetCol(iv)); el; ++el) {
            if (el.value() > 0) {
                if (model.GetRhs(el.index()).le != -kInf) n_down_locks++;
                if (model.GetRhs(el.index()).ri != kInf) n_up_locks++;
            }
            if (el.value() < 0) {
                if (model.GetRhs(el.index()).le != -kInf) n_up_locks++;
                if (model.GetRhs(el.index()).ri != kInf) n_down_locks++;
            }
        }

        if (n_up_locks == 0) {
            if (n_down_locks == 0) {
                // Todo it's either free or does not appear in constraints
                std::cerr << "Handle me\n";
                exit(0);
            } else {
                sol.x[iv] = Ceil(sol.x[iv]);
            }
        } else {
            if (n_down_locks == 0) {
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
