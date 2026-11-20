#pragma once

#include <array>

#include "reshala/heuristics/manager.h"
#include "reshala/milp/history.h"
#include "reshala/milp/mip_tracker.h"
#include "reshala/model/milp_model.h"

namespace reshala {

const Scalar kFsbMu = 1. / 6.;

struct BranchStats {
    Index n_lps = 0;
    Index n_lp_iters = 0;
};
std::ostream& operator<<(std::ostream& os, const BranchStats& stats);

class Branching {
   public:
    Branching(Reshala& ctx, History& hist)
        : ctx_(ctx),
          model_(ctx.GetModel()),
          mip_tracker_(ctx.GetMip()),
          heur_manager_(ctx.GetHeurMng()),
          hist_(hist) {}

    Index Branch(Node& parent, DualSimplex& ds);

    inline const Node& GetChild(Index i) const { return children_[i]; }

    Index FindBestChild() {
        return (children_[1].sol.status != LpStatus::kOptimal)
                   ? 0
                   : ((children_[0].sol.status != LpStatus::kOptimal)
                          ? 1
                          : ((children_[0].sol.y < children_[1].sol.y) ? 0 : 1));
    }

    inline const BranchStats& GetStats() const { return stats; }

   private:
    Reshala& ctx_;
    const MilpModel& model_;
    MipTracker& mip_tracker_;
    HeuristicManager& heur_manager_;
    History& hist_;

    BranchStats stats;

    std::array<Node, 2> children_;
    Index best_child_;
};

}  // namespace reshala
