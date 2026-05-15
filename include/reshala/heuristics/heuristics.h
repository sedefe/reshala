#pragma once

#include "reshala/heuristics/diving.h"

namespace reshala {

class Heuristics {
   public:
    Heuristics(MilpModel& model, MipState& mip_state);

    Solution Run(const Solution& relaxed);

   private:
    std::vector<std::unique_ptr<AbstractHeuristic>> heuristics;
    MilpModel& model_;
    MipState& mip_state_;

    Solution best_sol_;
};

}  // namespace reshala
