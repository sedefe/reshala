#pragma once

#include <memory>

#include "reshala/logging.h"
#include "reshala/lp/dual_simplex.h"
#include "reshala/milp/branching.h"
#include "reshala/milp/utils.h"

namespace reshala {

struct BnbStats {
    Index n_nodes = 0;
    Index n_dropped = 0;
};
std::ostream& operator<<(std::ostream& os, const BnbStats& stats);

class BnbSolver {
   public:
    BnbSolver(MilpModel& model, DualSimplex& ds, MipState& mip_state);

    void Solve(const Solution& relaxed);

    inline const BnbStats& GetStats() const { return stats; }

   private:
    MilpModel& model_;
    DualSimplex& ds_;
    MipState& mip_state_;

    BnbStats stats;

    std::vector<Node> nodes;

    Node curr_node;

    std::unique_ptr<AbstractBranching> branching_;

    void SolveRoot(Node& root);
    void UpdDual();

    void DebugPrint();
};

}  // namespace reshala
