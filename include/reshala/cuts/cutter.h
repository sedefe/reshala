#pragma once

#include "reshala/cuts/cut.h"
#include "reshala/cuts/generators/probing.h"
#include "reshala/milp/utils.h"

namespace reshala {

class Cutter {
   public:
    Cutter(MilpModel& model, Presolver& presolver, DualSimplex& ds, MipState& mip_state);

    void Run(Solution& relaxed);

   private:
    MilpModel& model_;
    DualSimplex& ds_;
    MipState& mip_state_;

    std::vector<std::unique_ptr<AbstractCg>> generators_;

    std::vector<Cut> cuts_;
    std::vector<Cut> fresh_cuts_;

    Index Generate(const Solution& relaxed);
    void CalcMetrics();
    Index Select();
    Index Add();
};

}  // namespace reshala
