#pragma once

#include "reshala/logging.h"
#include "reshala/lp/dual_simplex.h"
#include "reshala/milp/bnb.h"
#include "reshala/milp/branching.h"
#include "reshala/milp/utils.h"

namespace reshala {

class BnbSolver {
   public:
    BnbSolver(MilpModel& model, MipState& mip_state)
        : model_(model), mip_state_(mip_state), branching(model) {}

    Solution Solve(const Solution& sol);

   private:
    MipState& mip_state_;
    MilpModel& model_;
    std::vector<Node> nodes;
    Index n_nodes_ = 0;

    MostInfeasible branching;

    void SolveRoot();
    void UpdDual();

    void DebugPrint();
};

}  // namespace reshala
