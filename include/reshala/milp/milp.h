#pragma once

#include "reshala/milp/bnb.h"
#include "reshala/presolve/presolve.h"

namespace reshala {

class MilpSolver {
   public:
    MilpSolver(MilpModel& model);

    Solution Solve();

    MilpModel& model;
    MipState mip_state;

    Presolver presolver;
    BnbSolver bnb;

   private:
};

}  // namespace reshala
