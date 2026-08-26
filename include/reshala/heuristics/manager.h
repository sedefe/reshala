#pragma once

#include "reshala/heuristics/diving.h"
#include "reshala/heuristics/rounding.h"

namespace reshala {

enum class HeuristicTrigger { kRoot, kCut, kFsb, kNode };

struct HeurFreq {
    std::unique_ptr<AbstractHeuristic> h;
    Index freq;
};

class HeuristicManager {
   public:
    HeuristicManager(MipTracker& mip_tracker);
    void Run(HeuristicTrigger trigger, const MilpModel& model, const Solution& relaxed);

    void PrintStats(std::ostream& os) const;

   private:
    std::vector<HeurFreq> heuristics_;
    Rounding rounding;
    MipTracker& mip_tracker_;

    Index n_nodes_ = 0;
};

}  // namespace reshala
