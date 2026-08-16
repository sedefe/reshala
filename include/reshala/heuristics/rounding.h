#pragma once

#include "reshala/heuristics/abstract.h"
#include "reshala/presolve/presolve.h"

namespace reshala {

class Rounding : public AbstractHeuristic {
   public:
    Rounding(HeuristicType type) : AbstractHeuristic("Rounding", type) {}

   protected:
    Solution InternalRun(const MilpModel& model, const Solution& relaxed,
                         const MipTracker& mip_tracker);

   private:
};

}  // namespace reshala
