#pragma once

#include "reshala/cuts/cutter.h"
#include "reshala/heuristics/heuristics.h"
#include "reshala/milp/bnb.h"
#include "reshala/presolve/presolve.h"

namespace reshala {

class MilpSolver {
   public:
    MilpSolver(MilpModel& model);

    Solution Solve();

    MilpModel& model;
    DualSimplex ds;
    MipState mip_state;

    Presolver presolver;
    Cutter cutter;
    Heuristics heuristics;
    BnbSolver bnb;

    void PrintStats(std::ostream& os) const;

   private:
};

}  // namespace reshala
