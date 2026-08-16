#include "reshala/heuristics/manager.h"

namespace reshala {

HeuristicManager::HeuristicManager(MipTracker& mip_tracker) : mip_tracker_(mip_tracker) {
    std::vector<std::unique_ptr<AbstractHeuristic>> heuristics_;
    heuristics_.push_back(std::make_unique<Rounding>(HeuristicType::kFast));
    heuristics_.push_back(std::make_unique<Diving>(HeuristicType::kSlow, FixingType::kInts));

    for (auto& heur : heuristics_) {
        HeuristicType type = heur->type;
        heur_map_[type].push_back(std::move(heur));
    }
}

void HeuristicManager::Run(HeuristicTrigger trigger, const MilpModel& model,
                           const Solution& relaxed) {
    if (relaxed.status != LpStatus::kOptimal) return;
    n_tries_++;

    if (trigger == HeuristicTrigger::kRoot or trigger == HeuristicTrigger::kCut) {
        for (auto& [type, heuristics] : heur_map_) {
            for (auto& h : heuristics) {
                Solution sol = h->Run(model, relaxed, mip_tracker_);
                if (mip_tracker_.TestPrimal(sol)) {
                    ReportNewPrimal(h->GetName(), sol.y);
                }
                if (mip_tracker_.Converged()) {
                    break;
                }
            }
        }
    } else if (trigger == HeuristicTrigger::kBnb) {
        for (auto& h : heur_map_[HeuristicType::kFast]) {
            Solution sol = h->Run(model, relaxed, mip_tracker_);
            if (mip_tracker_.TestPrimal(sol)) {
                ReportNewPrimal(h->GetName(), sol.y);
            }
            if (mip_tracker_.Converged()) {
                break;
            }
        }
    }
}

}  // namespace reshala
