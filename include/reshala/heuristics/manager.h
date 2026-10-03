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
    HeuristicManager(Reshala& ctx);
    void RunStartHeus(HeuristicTrigger trigger, const Solution& relaxed);

    void PrintStats(std::ostream& os) const;

   private:
    Reshala& ctx_;
    std::vector<HeurFreq> start_heus_;
    Rounding rounding_;

    Index n_nodes_ = 0;
};

}  // namespace reshala
