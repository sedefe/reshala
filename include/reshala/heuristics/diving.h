#pragma once

#include "reshala/heuristics/abstract.h"
#include "reshala/presolve/presolve.h"

namespace reshala {

class Diving : public AbstractHeuristic {
   public:
    Diving(Reshala& ctx, FixingType fixing_type)
        : AbstractHeuristic(ctx, "Diving-" + FixingType2Str(fixing_type)),
          fixing_type_(fixing_type) {}

   protected:
    Solution InternalRun(const Solution& relaxed);
    FixingType fixing_type_;

   private:
    Index GetCandidate(const MilpModel& model, const Solution& relaxed, const Solution& sol);
};

}  // namespace reshala
