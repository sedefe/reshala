#pragma once

#include "reshala/cuts/cutter.h"
#include "reshala/heuristics/manager.h"
#include "reshala/milp/bnb.h"
#include "reshala/presolve/presolve.h"

namespace reshala {

class MilpSolver {
   public:
    MilpSolver(MilpModel& model);

    Solution Solve();

    MilpModel& model;
    DualSimplex ds;
    MipTracker mip_tracker;

    Presolver presolver;
    HeuristicManager heur_manager;
    Cutter cutter;
    BnbSolver bnb;

    void PrintStats(std::ostream& os) const;

   private:
};

}  // namespace reshala
