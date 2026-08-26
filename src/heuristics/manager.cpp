#include "reshala/heuristics/manager.h"

namespace reshala {

HeuristicManager::HeuristicManager(MipTracker& mip_tracker) : mip_tracker_(mip_tracker) {
    heuristics_.push_back({std::make_unique<Rounding>(), 1});
    heuristics_.push_back({std::make_unique<Diving>(FixingType::kAll), 10});
    heuristics_.push_back({std::make_unique<Diving>(FixingType::kInts), 15});
    heuristics_.push_back({std::make_unique<Diving>(FixingType::kNone), 50});
}

void HeuristicManager::Run(HeuristicTrigger trigger, const MilpModel& model,
                           const Solution& relaxed) {
    if (relaxed.status != LpStatus::kOptimal) return;
    if (relaxed.y >= mip_tracker_.GetCutoff()) return;

    switch (trigger) {
        case HeuristicTrigger::kRoot:
        case HeuristicTrigger::kCut:
            for (auto& [h, freq] : heuristics_) {
                h->Run(model, relaxed, mip_tracker_);
                if (mip_tracker_.Converged()) {
                    break;
                }
            }
            break;
        case HeuristicTrigger::kNode:
            for (auto& [h, freq] : heuristics_) {
                if (n_nodes_ % freq == 0) {
                    h->Run(model, relaxed, mip_tracker_);
                    if (mip_tracker_.Converged()) {
                        break;
                    }
                }
            }
            n_nodes_++;
            break;
        case HeuristicTrigger::kFsb:
            rounding.Run(model, relaxed, mip_tracker_);
            if (mip_tracker_.Converged()) {
                break;
            }
            break;
        default:
            assert(false && "Unknown heuristic trigger");
    }
}

void HeuristicManager::PrintStats(std::ostream& os) const {
    os << "Heuristics:\n";
    for (auto& [h, freq] : heuristics_) {
        os << "\t" << std::setw(12) << h->GetName() << ": " << h->stats << "\n";
    }
}

}  // namespace reshala
