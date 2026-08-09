#include "reshala/heuristics/heuristics.h"

namespace reshala {

Heuristics::Heuristics(const MilpModel& model, MipTracker& mip_tracker)
    : model_(model), mip_tracker_(mip_tracker) {
    heuristics_.push_back(std::make_unique<Rounding>());
    heuristics_.push_back(std::make_unique<Diving>(FixingType::kInts));
}

void Heuristics::Run(const Solution& relaxed) {
    for (auto& h : heuristics_) {
        Solution sol = h->Run(model_, relaxed, mip_tracker_);
        mip_tracker_.TestPrimal(sol);
        if (best_sol_.y > sol.y) {
            best_sol_ = sol;
        }

        if (mip_tracker_.Converged()) {
            break;
        }
    }
}

}  // namespace reshala
