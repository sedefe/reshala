#pragma once

#include "reshala/heuristics/abstract.h"
#include "reshala/presolve/presolve.h"

namespace reshala {

class Diving : public AbstractHeuristic {
   public:
    Diving(RoundingType type)
        : AbstractHeuristic("Diving-" + RoundingType2Str(type)), type_(type) {}

   protected:
    Solution InternalRun(const MilpModel& model, const Solution& relaxed, const MipState& mip_state);
    RoundingType type_;

   private:
    Index GetCandidate(const MilpModel& model, const Solution& relaxed, const Solution& sol);
};

}  // namespace reshala
