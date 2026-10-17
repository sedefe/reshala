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

    LpStatus Presolve(bool verbose, bool hard = true);
    Solution Postsolve(const Solution&);

    inline const ModelTracker& GetTracker() const { return *tracker_; }

   private:
    MilpModel& model_;
    std::unique_ptr<ModelTracker> tracker_;
    RuleResult status_;

    std::vector<std::unique_ptr<RowRule>> row_rules_;
    std::vector<std::unique_ptr<ColRule>> col_rules_;
    std::vector<std::unique_ptr<HardRule>> hard_rules_;

    inline bool Solved() {
        return (status_ == RuleResult::kInfeasible) or (model_.GetNVars() == 0);
    }

    RuleResult SimplePresolve(bool verbose);

    void PrintHeader() const;
    void PrintStat(const Rule&, const PresolveStat& stat) const;
};

}  // namespace reshala
