#pragma once

#include "reshala/cuts/cutter.h"
#include "reshala/heuristics/manager.h"
#include "reshala/milp/bnb.h"
#include "reshala/presolve/presolve.h"
#include "reshala/reshala.h"

namespace reshala {

class MilpSolver {
   public:
    MilpSolver(Reshala& ctx);

    Solution Solve();

   private:
    Reshala& ctx_;
};

}  // namespace reshala
