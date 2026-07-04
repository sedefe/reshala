#pragma once

#include "reshala/cuts/cut.h"
#include "reshala/lp/dual_simplex.h"

namespace reshala {

class AbstractCg {
   public:
    AbstractCg(const std::string& name, MilpModel& model, DualSimplex& ds)
        : name_(name), model_(model), ds_(ds) {}
    virtual ~AbstractCg() = default;

    const std::string GetName() const { return name_; }
    virtual Index Generate(const Solution& sol, std::vector<Cut>& dst) = 0;

   protected:
    const std::string name_;
    MilpModel& model_;
    DualSimplex& ds_;
};

}  // namespace reshala
