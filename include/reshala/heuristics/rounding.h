#pragma once

#include "reshala/heuristics/abstract.h"
#include "reshala/presolve/presolve.h"

namespace reshala {

class Rounding : public AbstractHeuristic {
   public:
    Rounding(Reshala& ctx) : AbstractHeuristic(ctx, "Rounding") {}

   protected:
    Solution InternalRun(const Solution& relaxed);

   private:
};

}  // namespace reshala
