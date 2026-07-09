#pragma once

#include "reshala/cuts/generators/abstract_cg.h"
#include "reshala/presolve/presolve.h"

namespace reshala {

class ProbingCg : public AbstractCg {
   public:
    ProbingCg(MilpModel& model, Presolver& presolver, DualSimplex& ds)
        : AbstractCg("Probing", model, ds), impls_(presolver.GetTracker().GetImplications()) {}

    void Generate(const Solution& sol, std::vector<Cut>& dst) override;

   private:
    const std::vector<Implication>& impls_;
};

}  // namespace reshala
