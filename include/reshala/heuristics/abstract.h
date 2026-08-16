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

    Solution Run(const MilpModel& model, const Solution& relaxed, const MipTracker& mip_tracker,
                 bool verbose) {
        if (verbose) std::cout << "Running " << name_ << "\n";

        auto [sol, t_heur] = MEASURE_TIME(InternalRun(model, relaxed, mip_tracker));
        if (verbose) {
            std::cout << "Finished in " << t_heur << " ms: ";
            if (sol.status == LpStatus::kOptimal) {
                std::cout << "found solution " << FMT(-10, 5) << sol.y << "\n";
            } else {
                std::cout << "did not find any solution\n";
            }
        }
        return sol;
    }

    HeuristicType type;

   protected:
    const std::string name_;
    virtual Solution InternalRun(const MilpModel& model, const Solution& relaxation,
                                 const MipTracker& mip_tracker) = 0;
};

}  // namespace reshala
