#pragma once

#include "reshala/model/milp_model.h"
#include "reshala/presolve/model_tracker.h"
#include "reshala/presolve/rules.h"
#include "reshala/reshala.h"

namespace reshala {

class Presolver {
   public:
    Presolver(Reshala& ctx);

    void SimplePresolve(bool verbose);
    LpStatus Presolve(MilpModel& model, bool verbose, bool hard = true);
    Solution Postsolve(const Solution&);

    inline const ModelTracker& GetTracker() const { return *tracker_; }

   private:
    Reshala& ctx_;

    std::unique_ptr<ModelTracker> tracker_;
    RuleResult status_;

    std::vector<std::unique_ptr<Rule>> simple_rules_;
    std::vector<std::unique_ptr<Rule>> hard_rules_;

    void PrintHeader() const;
    void PrintStat(const Rule&, const PresolveStat& stat) const;
};

}  // namespace reshala
