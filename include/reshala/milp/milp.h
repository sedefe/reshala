#pragma once

#include "reshala/milp/bnb.h"
#include "reshala/presolve/presolve.h"

namespace reshala {

class MilpSolver {
   public:
    MilpSolver(MilpModel& model);

    Solution Solve();

    MipState mip_state;
    MilpModel& model;

    Presolver presolver;
    DualSimplex ds;
    BnbSolver bnb;

   private:
};

}  // namespace reshala
