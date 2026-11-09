#include "reshala/reshala.h"

#include "reshala/cuts/cutter.h"
#include "reshala/io/io.h"
#include "reshala/milp/milp.h"

namespace reshala {

Reshala::Reshala()
    : io_(std::make_unique<Io>(*this)),
      presolver_(std::make_unique<Presolver>(*this)),
      ds_(std::make_unique<DualSimplex>(*this)),
      mip_tracker_(std::make_unique<MipTracker>(*this)),
      heur_manager_(std::make_unique<HeuristicManager>(*this)),
      cutter_(std::make_unique<Cutter>(*this)),
      bnb_(std::make_unique<BnbSolver>(*this)),
      milp_(std::make_unique<MilpSolver>(*this)) {}

Reshala::~Reshala() = default;

FileReadStatus Reshala::Read(const char* path) {
    auto res = io_->Read(path);
    return res;
}

Solution Reshala::Solve() { return milp_->Solve(); }

void Reshala::PrintModel(std::ostream& os, const MilpModel& model) const {
    io_->PrintModel(os, model);
}

void Reshala::PrintSolution(std::ostream& os, const std::vector<Scalar>& x) const {
    const auto& names = io_->GetNames();
    assert(names.vars.Size() == x.size());
    for (Index iv = 0; iv < x.size(); iv++) {
        if (!IsZero(x[iv])) {
            os << names.vars.GetName(iv) << ": " << x[iv] << " ";
        }
    }
    os << "\n";
}

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
