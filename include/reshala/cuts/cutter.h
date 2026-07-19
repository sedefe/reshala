#pragma once

#include "reshala/cuts/config.h"
#include "reshala/cuts/cut.h"
#include "reshala/cuts/generators/probing.h"
#include "reshala/milp/utils.h"

namespace reshala {

class Cutter {
   public:
    Cutter(MilpModel& model, const Presolver& presolver, DualSimplex& ds, MipState& mip_state);

    void Run(Solution& sol);

   private:
    MilpModel& model_;
    const Presolver& presolver_;
    DualSimplex& ds_;
    MipState& mip_state_;
    Solution sol_;

    Index n_round_ = 0;

    Index max_cuts_;
    Index max_support_;

    std::vector<Cut> pool_;

    Index Generate(const Solution& sol);
    void Estimate();
    void Filter();
    Index Select();
    Index Add();

    bool AnyNewSelected();
};

}  // namespace reshala
