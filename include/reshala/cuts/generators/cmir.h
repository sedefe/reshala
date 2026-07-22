#pragma once

#include "reshala/cuts/generators/abstract_cg.h"
#include "reshala/presolve/presolve.h"

namespace reshala {

class CmirCg : public AbstractCg {
   public:
    CmirCg(MilpModel& model, DualSimplex& ds) : AbstractCg("Cmir", model, ds) {}

    void Generate(const Solution& sol, std::vector<Cut>& dst) override;

   private:
};

}  // namespace reshala
