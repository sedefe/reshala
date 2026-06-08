#pragma once

#include "reshala/milp/utils.h"
#include "reshala/model/milp_model.h"

namespace reshala {

const Scalar kFsbMu = 1. / 6.;

class AbstractBranching {
   public:
    AbstractBranching(MilpModel& model) : model_(model) {}
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
    std::array<Node, 2> children_;
    Index best_child_;
};

class MostInfeasible : public AbstractBranching {
   public:
    MostInfeasible(MilpModel& model) : AbstractBranching(model) {}
    Index Branch(Node& parent, DualSimplex& ds) override;
};

class FullStrong : public AbstractBranching {
   public:
    FullStrong(MilpModel& model) : AbstractBranching(model) {}
    Index Branch(Node& parent, DualSimplex& ds) override;
};

class FullStrongDomProp : public AbstractBranching {
   public:
    FullStrongDomProp(MilpModel& model) : AbstractBranching(model) {}
    Index Branch(Node& parent, DualSimplex& ds) override;
};

}  // namespace reshala
