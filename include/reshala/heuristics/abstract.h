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
    AbstractHeuristic(Reshala& ctx, const std::string& name) : ctx_(ctx), name_(name) {}
    virtual ~AbstractHeuristic() = default;
    const std::string& GetName() const { return name_; }

    void Run(const Solution& relaxed) {
        stats.n_called++;

        auto [sol, t_heur] = MEASURE_TIME(InternalRun(relaxed));
        stats.time += t_heur;

        if (sol.status == LpStatus::kOptimal) {
            stats.n_found++;
            if (ctx_.GetMip().TestPrimal(sol)) {
                ReportNewPrimal(GetName(), sol.y);
                stats.n_improved++;
            }
        }
    }

    HeurStats stats;

   protected:
    Reshala& ctx_;
    const std::string name_;
    virtual Solution InternalRun(const Solution& relaxation) = 0;
};

}  // namespace reshala
