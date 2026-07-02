#pragma once

#include <array>

#include "reshala/milp/utils.h"
#include "reshala/model/milp_model.h"

namespace reshala {

const Scalar kFsbMu = 1. / 6.;

class AbstractBranching {
   public:
    AbstractBranching(MilpModel& model, MipState& mip_state)
        : model_(model), mip_state_(mip_state) {}
    virtual Index Branch(Node& parent, DualSimplex& ds) = 0;

    inline const Node& GetChild(Index i) const { return children_[i]; }

    Index FindBestChild() {
        return (children_[1].sol.status != LpStatus::kOptimal)
                   ? 0
                   : ((children_[0].sol.status != LpStatus::kOptimal)
                          ? 1
                          : ((children_[0].sol.y < children_[1].sol.y) ? 0 : 1));
    }

   protected:
    MilpModel& model_;
    MipState& mip_state_;
    std::array<Node, 2> children_;
    Index best_child_;
};

class MostInfeasible : public AbstractBranching {
   public:
    MostInfeasible(MilpModel& model, MipState& mip_state) : AbstractBranching(model, mip_state) {}
    Index Branch(Node& parent, DualSimplex& ds) override;
};

class FullStrong : public AbstractBranching {
   public:
    FullStrong(MilpModel& model, MipState& mip_state) : AbstractBranching(model, mip_state) {}
    Index Branch(Node& parent, DualSimplex& ds) override;
};

}  // namespace reshala
