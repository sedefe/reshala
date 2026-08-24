#include "reshala/heuristics/manager.h"

namespace reshala {

HeuristicManager::HeuristicManager(MipTracker& mip_tracker) : mip_tracker_(mip_tracker) {
    std::vector<std::unique_ptr<AbstractHeuristic>> heuristics_;
    heuristics_.push_back(std::make_unique<Rounding>(HeuristicType::kFast));
    heuristics_.push_back(std::make_unique<Diving>(HeuristicType::kSlow, FixingType::kAll));
    heuristics_.push_back(std::make_unique<Diving>(HeuristicType::kSlow, FixingType::kInts));
    heuristics_.push_back(std::make_unique<Diving>(HeuristicType::kSlow, FixingType::kNone));

    for (auto& heur : heuristics_) {
        HeuristicType type = heur->type;
        heur_map_[type].push_back(std::move(heur));
    }
}

void HeuristicManager::Run(HeuristicTrigger trigger, const MilpModel& model,
                           const Solution& relaxed) {
    if (relaxed.status != LpStatus::kOptimal) return;
    if (relaxed.y >= mip_tracker_.GetCutoff()) return;
    n_tries_++;

    switch (trigger) {
        case HeuristicTrigger::kRoot:
        case HeuristicTrigger::kCut:
            for (auto& [type, heuristics] : heur_map_) {
                for (auto& h : heuristics) {
                    h->Run(model, relaxed, mip_tracker_);
                    if (mip_tracker_.Converged()) {
                        break;
                    }
                }
            }
            break;
        case HeuristicTrigger::kFsb:
        case HeuristicTrigger::kNode:
            for (auto& h : heur_map_[HeuristicType::kFast]) {
                h->Run(model, relaxed, mip_tracker_);
                if (mip_tracker_.Converged()) {
                    break;
                }
            }
            break;
        default:
            assert(false && "Unknown heuristic trigger");
    }
}

void HeuristicManager::PrintStats(std::ostream& os) const {
    os << "Heuristics:\n";
    for (auto& [type, heuristics] : heur_map_) {
        for (auto& h : heuristics) {
            os << "\t" << std::setw(12) << h->GetName() << ": " << h->stats << "\n";
        }
    }
}

}  // namespace reshala
