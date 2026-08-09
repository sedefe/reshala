#pragma once

#include "reshala/heuristics/abstract.h"
#include "reshala/presolve/presolve.h"

namespace reshala {

class Diving : public AbstractHeuristic {
   public:
    Diving(FixingType type)
        : AbstractHeuristic("Diving-" + FixingType2Str(type)), type_(type) {}

   protected:
    Solution InternalRun(const MilpModel& model, const Solution& relaxed, const MipTracker& mip_tracker);
    FixingType type_;

   private:
    Index GetCandidate(const MilpModel& model, const Solution& relaxed, const Solution& sol);
};

}  // namespace reshala
