#pragma once

#include "reshala/heuristics/utils.h"
#include "reshala/milp/mip_tracker.h"

namespace reshala {

struct HeurStats {
    Index n_called = 0;
    Index n_found = 0;
    Index n_improved = 0;
    Scalar time = 0.;
};
inline std::ostream& operator<<(std::ostream& os, const HeurStats& stats) {
    os << stats.n_called << " calls (" << FMT(0, 3) << stats.time << " ms), found " << stats.n_found
       << ", improved " << stats.n_improved;
    return os;
}

class AbstractHeuristic {
   public:
    AbstractHeuristic(const std::string& name) : name_(name) {}
    virtual ~AbstractHeuristic() = default;
    const std::string& GetName() const { return name_; }

    void Run(const MilpModel& model, const Solution& relaxed, MipTracker& mip_tracker) {
        stats.n_called++;

        auto [sol, t_heur] = MEASURE_TIME(InternalRun(model, relaxed, mip_tracker));
        stats.time += t_heur;

        if (sol.status == LpStatus::kOptimal) {
            stats.n_found++;
            if (mip_tracker.TestPrimal(sol)) {
                ReportNewPrimal(GetName(), sol.y);
                stats.n_improved++;
            }
        }
    }

    HeurStats stats;

   protected:
    const std::string name_;
    virtual Solution InternalRun(const MilpModel& model, const Solution& relaxation,
                                 const MipTracker& mip_tracker) = 0;
};

}  // namespace reshala
