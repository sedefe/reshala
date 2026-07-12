#include "reshala/heuristics/heuristics.h"

namespace reshala {

Heuristics::Heuristics(const MilpModel& model, MipState& mip_state)
    : model_(model), mip_state_(mip_state) {
    heuristics_.push_back(std::make_unique<Diving>(RoundingType::kInts));
}

void Heuristics::Run(const Solution& relaxed) {
    for (auto& h : heuristics_) {
        Solution sol = h->Run(model_, relaxed, mip_state_);
        mip_state_.TestPrimal(sol);
        if (best_sol_.y > sol.y) {
            best_sol_ = sol;
        }

        if (mip_state_.Converged()) {
            break;
        }
    }
}

}  // namespace reshala
