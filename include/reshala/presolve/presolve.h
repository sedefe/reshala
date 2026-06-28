#pragma once

#include <map>
#include <memory>

#include "reshala/model/milp_model.h"
#include "reshala/presolve/model_tracker.h"
#include "reshala/presolve/rules.h"

namespace reshala {

class Presolver {
   public:
    Presolver(MilpModel& model);

    LpStatus Presolve();
    Solution Postsolve(const Solution&);

    inline const ModelTracker& GetTracker() { return tracker_; }

   private:
    MilpModel& model_;
    ModelTracker tracker_;

    std::map<RuleType, std::vector<std::unique_ptr<Rule>>> rule_map_;

    void PrintHeader() const;
    void PrintStat(const Rule&, const PresolveStat& stat) const;
};

}  // namespace reshala
