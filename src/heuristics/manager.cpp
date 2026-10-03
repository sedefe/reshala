#include "reshala/heuristics/manager.h"

namespace reshala {

HeuristicManager::HeuristicManager(Reshala& ctx) : ctx_(ctx), rounding_(ctx) {
    start_heus_.push_back({std::make_unique<Rounding>(ctx), 1});
    start_heus_.push_back({std::make_unique<Diving>(ctx, FixingType::kAll), 10});
    start_heus_.push_back({std::make_unique<Diving>(ctx, FixingType::kInts), 15});
    start_heus_.push_back({std::make_unique<Diving>(ctx, FixingType::kNone), 50});
}

void HeuristicManager::RunStartHeus(HeuristicTrigger trigger, const Solution& relaxed) {
    if (relaxed.status != LpStatus::kOptimal) return;
    if (relaxed.y >= ctx_.GetMip().GetCutoff()) return;

    MipTracker& mip_tracker = ctx_.GetMip();

    switch (trigger) {
        case HeuristicTrigger::kRoot:
        case HeuristicTrigger::kCut:
            for (auto& [h, freq] : start_heus_) {
                h->Run(relaxed);
                if (mip_tracker.Converged()) {
                    break;
                }
            }
            break;
        case HeuristicTrigger::kNode:
            for (auto& [h, freq] : start_heus_) {
                if (n_nodes_ % freq == 0) {
                    h->Run(relaxed);
                    if (mip_tracker.Converged()) {
                        break;
                    }
                }
            }
            n_nodes_++;
            break;
        case HeuristicTrigger::kFsb:
            rounding_.Run(relaxed);
            if (mip_tracker.Converged()) {
                break;
            }
            break;
        default:
            assert(false && "Unknown heuristic trigger");
    }
}

void HeuristicManager::PrintStats(std::ostream& os) const {
    os << "Heuristics:\n";
    for (auto& [h, freq] : start_heus_) {
        os << "\t" << std::setw(12) << h->GetName() << ": " << h->stats << "\n";
    }
}

}  // namespace reshala
