#pragma once

#include "reshala/heuristics/utils.h"
#include "reshala/milp/mip_tracker.h"

namespace reshala {

enum class HeuristicType { kFast, kSlow };

class AbstractHeuristic {
   public:
    AbstractHeuristic(const std::string& name, HeuristicType t) : name_(name), type(t) {}
    virtual ~AbstractHeuristic() = default;
    const std::string& GetName() const { return name_; }

    Solution Run(const MilpModel& model, const Solution& relaxed, const MipTracker& mip_tracker) {
        auto [sol, t_heur] = MEASURE_TIME(InternalRun(model, relaxed, mip_tracker));
        // Todo keep time stats
        return sol;
    }

    HeuristicType type;

   protected:
    const std::string name_;
    virtual Solution InternalRun(const MilpModel& model, const Solution& relaxation,
                                 const MipTracker& mip_tracker) = 0;
};

}  // namespace reshala
