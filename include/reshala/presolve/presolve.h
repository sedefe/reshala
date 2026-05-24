#pragma once

#include <memory>

#include "reshala/model/milp_model.h"
#include "reshala/presolve/model_tracker.h"
#include "reshala/presolve/rules.h"
#include "reshala/presolve/transforms.h"

namespace reshala {

class Presolver {
   public:
    Presolver(MilpModel& model);

    LpStatus Presolve();
    Solution Postsolve(const Solution&);

   private:
    MilpModel& model_;
    ModelTracker tracker_;

    std::vector<std::unique_ptr<Rule>> rules_;
    std::vector<std::unique_ptr<Transform>> transforms_;

    void PrintHeader() const;
    void PrintStat(const Rule&, const PresolveStat& stat) const;
};

}  // namespace reshala
