#include "reshala/heuristics/heuristics.h"

namespace reshala {

Heuristics::Heuristics(MilpModel& model, MipState& mip_state)
    : model_(model), mip_state_(mip_state) {
    heuristics.push_back(std::make_unique<Diving>(RoundingType::kInts));
}

Solution Heuristics::Run(const Solution& relaxed) {
    for (auto& h : heuristics) {
        MilpModel model_copy = model_;

        Solution sol = h->Run(model_copy, relaxed, mip_state_);
        if (best_sol_.y < sol.y) {
            best_sol_ = sol;
        }

        if (mip_state_.Converged()) {
            break;
        }
    }

    return best_sol_;
}

}  // namespace reshala
