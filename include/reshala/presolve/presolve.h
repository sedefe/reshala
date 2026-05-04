#pragma once

#include "reshala/model/milp_model.h"

namespace reshala {

class Presolver {
   public:
    Presolver(MilpModel& model);

    void Presolve();

   private:
    MilpModel& model_;
};

}  // namespace reshala
