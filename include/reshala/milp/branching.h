#pragma once

#include "reshala/model/milp_model.h"

namespace reshala {

class AbstractBranching {
   public:
    AbstractBranching(const MilpModel& model) : model_(model) {}
    virtual Index Branch(const Solution& sol) = 0;

   protected:
    const MilpModel& model_;
};

class MostInfeasible : public AbstractBranching {
   public:
    MostInfeasible(const MilpModel& model) : AbstractBranching(model) {}
    Index Branch(const Solution& sol) override;
};

class FullStrong : public AbstractBranching {
   public:
    FullStrong(const MilpModel& model) : AbstractBranching(model) {}
    Index Branch(const Solution& sol) override;
};

}  // namespace reshala
