#pragma once

#include "reshala/heuristics/diving.h"
#include "reshala/heuristics/rounding.h"

namespace reshala {

enum class HeuristicTrigger { kRoot, kCut, kFsb, kNode };

class HeuristicManager {
   public:
    HeuristicManager(MipTracker& mip_tracker);
    void Run(HeuristicTrigger trigger, const MilpModel& model, const Solution& relaxed);

    void PrintStats(std::ostream& os) const;

   private:
    std::map<HeuristicType, std::vector<std::unique_ptr<AbstractHeuristic>>> heur_map_;
    MipTracker& mip_tracker_;

    Index n_tries_ = 0;
};

}  // namespace reshala
