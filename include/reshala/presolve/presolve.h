#pragma once

#include "reshala/model/milp_model.h"

namespace reshala {

class Presolver {
   public:
    Presolver(MilpModel& model_);

   private:
    MilpModel& model;
};

}  // namespace reshala
