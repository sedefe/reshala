#pragma once

#include "reshala/heuristics/diving.h"

namespace reshala {

class Heuristics {
   public:
    Heuristics(const MilpModel& model, MipTracker& mip_tracker);

    void Run(const Solution& relaxed);

   private:
    std::vector<std::unique_ptr<AbstractHeuristic>> heuristics_;
    const MilpModel& model_;
    MipTracker& mip_tracker_;

    Solution best_sol_;
};

}  // namespace reshala
