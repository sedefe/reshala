#pragma once

#include "reshala/model/solution.h"

namespace reshala {

class AbstractBranching {
   public:
    virtual Index Branch(const Solution& sol) = 0;
};

class MostInfeasible : public AbstractBranching {
   public:
    Index Branch(const Solution& sol) override;
};

class FullStrong : public AbstractBranching {
   public:
    Index Branch(const Solution& sol) override;
};

}  // namespace reshala
