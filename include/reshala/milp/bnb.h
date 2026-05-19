#pragma once

#include "reshala/logging.h"
#include "reshala/lp/dual_simplex.h"
#include "reshala/milp/bnb.h"
#include "reshala/milp/branching.h"
#include "reshala/milp/utils.h"

namespace reshala {

class BnbSolver {
   public:
    BnbSolver(MilpModel& model, DualSimplex& ds, MipState& mip_state)
        : model_(model), ds_(ds), mip_state_(mip_state), branching(model) {}

    void Solve(const Solution& relaxed);

   private:
    MilpModel& model_;
    DualSimplex& ds_;
    MipState& mip_state_;
    std::vector<Node> nodes;
    Index n_nodes_ = 0;

    Node curr_node;

    MostInfeasible branching;

    void SolveRoot();
    void UpdDual();

    void DebugPrint();
};

}  // namespace reshala
