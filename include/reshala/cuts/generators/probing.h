#pragma once

#include "reshala/cuts/generators/abstract_cg.h"
#include "reshala/presolve/presolve.h"

namespace reshala {

class ProbingCg : public AbstractCg {
   public:
    ProbingCg(Reshala& ctx)
        : AbstractCg(ctx, "Probing"), impls_(ctx.GetPresolver().GetTracker().GetImplications()) {}

    void Generate(const Solution& sol, std::vector<Cut>& dst) override;

   private:
    const std::vector<Implication>& impls_;
};

}  // namespace reshala
