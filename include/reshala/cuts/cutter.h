#pragma once

#include "reshala/cuts/cut.h"
#include "reshala/cuts/generators/probing.h"
#include "reshala/milp/utils.h"

namespace reshala {

class Cutter {
   public:
    Cutter(MilpModel& model, const Presolver& presolver, DualSimplex& ds, MipState& mip_state);

    void Run(Solution& sol);

   private:
    const Index kMaxRounds = 5;
    const Index kMaxAge = 3;
    const Scalar kMaxRelSupport = 0.01;
    const Scalar kMinRelQuality = 0.1;
    const Scalar kMaxCutsFactor = 1.0;
    const Scalar kThdCos2 = 1.0 / 100;

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
