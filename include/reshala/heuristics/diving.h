#pragma once

#include "reshala/heuristics/abstract.h"
#include "reshala/presolve/presolve.h"

namespace reshala {

class Diving : public AbstractHeuristic {
   public:
    Diving(HeuristicType type, FixingType fixing_type)
        : AbstractHeuristic("Diving-" + FixingType2Str(fixing_type), type),
          fixing_type_(fixing_type) {}

   protected:
    Solution InternalRun(const MilpModel& model, const Solution& relaxed,
                         const MipTracker& mip_tracker);
    FixingType fixing_type_;

   private:
    Index GetCandidate(const MilpModel& model, const Solution& relaxed, const Solution& sol);
};

}  // namespace reshala
