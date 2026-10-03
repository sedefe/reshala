#pragma once

#include "reshala/cuts/config.h"
#include "reshala/cuts/cut.h"
#include "reshala/lp/dual_simplex.h"

namespace reshala {

class AbstractCg {
   public:
    AbstractCg(Reshala& ctx, const std::string& name) : ctx_(ctx), name_(name) {}
    virtual ~AbstractCg() = default;

    const std::string GetName() const { return name_; }
    virtual void Generate(const Solution& sol, std::vector<Cut>& dst) = 0;

   protected:
    Reshala& ctx_;
    const std::string name_;
};

}  // namespace reshala
