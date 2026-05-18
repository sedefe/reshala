#pragma once

#include "reshala/heuristics/utils.h"
#include "reshala/milp/utils.h"

namespace reshala {

class AbstractHeuristic {
   public:
    AbstractHeuristic(const std::string& name) : name_(name) {}
    virtual ~AbstractHeuristic() = default;
    const std::string& GetName() const { return name_; }

    Solution Run(MilpModel& model, const Solution& relaxed, const MipState& mip_state) {
        std::cout << "Running " << name_ << "\n";

        auto [sol, duration] = MEASURE_TIME(InternalRun(model, relaxed, mip_state));
        std::cout << "Finished in " << duration.count() / 1e3 << " ms\n";

        if (sol.status == LpStatus::kOptimal) {
            std::cout << "Found solution " << FMT(-10, 5) << sol.y << "\n";
        } else {
            std::cout << "Did not find any solution\n";
        }
        return sol;
    }

   protected:
    const std::string name_;
    virtual Solution InternalRun(MilpModel& model, const Solution& relaxation,
                                 const MipState& mip_state) = 0;
};

}  // namespace reshala
