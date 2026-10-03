#include "reshala/reshala.h"

#include "reshala/cuts/cutter.h"
#include "reshala/milp/milp.h"

namespace reshala {

Reshala::Reshala()
    : io_(model_),
      presolver_(std::make_unique<Presolver>(model_)),
      ds_(std::make_unique<DualSimplex>(*this)),
      mip_tracker_(std::make_unique<MipTracker>(*this)),
      heur_manager_(std::make_unique<HeuristicManager>(*this)),
      cutter_(std::make_unique<Cutter>(*this)),
      bnb_(std::make_unique<BnbSolver>(*this)),
      milp_(std::make_unique<MilpSolver>(*this)) {}

Reshala::~Reshala() = default;

FileReadStatus Reshala::Read(const char* path) {
    auto res = io_.Read(path);
    return res;
}

Solution Reshala::Solve() { return milp_->Solve(); }

void Reshala::PrintStats(std::ostream& os) const {
    os << "=== Stats ===\n";
    os << ds_->GetLina().GetStats();
    os << ds_->GetStats();
    os << ds_->GetScaling().stats;
    heur_manager_->PrintStats(os);
    os << cutter_->GetStats();
    os << bnb_->GetStats();
}

}  // namespace reshala
