#pragma once

#include <optional>

#include "reshala/heuristics/heuristics.h"
#include "reshala/milp/bnb.h"
#include "reshala/presolve/presolve.h"

namespace reshala {

class MilpSolver {
   public:
    MilpSolver(MilpModel& model);

    Solution Solve();

    MilpModel& model;
    std::optional<DualSimplex> ds;  // Will be constructed after presolving
    MipState mip_state;

    Presolver presolver;
    Heuristics heuristics;
    BnbSolver bnb;

   private:
};

}  // namespace reshala
