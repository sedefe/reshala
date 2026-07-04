#include "reshala/cuts/cutter.h"

namespace reshala {

Cutter::Cutter(MilpModel& model, Presolver& presolver, DualSimplex& ds, MipState& mip_state)
    : model_(model), ds_(ds), mip_state_(mip_state) {
    generators_.push_back(std::make_unique<ProbingCg>(model, presolver, ds));
}

void Cutter::Run(Solution& relaxed) {
    auto n_generated = Generate(relaxed);
    CalcMetrics();
    auto n_selected = Select();
    auto n_added = Add();
    std::cout << "Added " << n_added << " cuts, size: " << model_.GetNCons() << " x "
              << model_.GetNVars() << "\n";

    if (n_added > 0) {
        ds_.SetModel(model_);
        relaxed = ds_.Solve(false);
        mip_state_.UpdDual(relaxed.y);
    }
    std::cout << "After cutter: " << relaxed.y << "\n";
}

Index Cutter::Generate(const Solution& relaxed) {
    fresh_cuts_.clear();
    for (auto& cg : generators_) {
        Index k = fresh_cuts_.size();
        cg->Generate(relaxed, fresh_cuts_);
        k = fresh_cuts_.size() - k;
        std::cout << cg->GetName() << " generated " << k << " cuts\n";
    }
    return fresh_cuts_.size();
}

void Cutter::CalcMetrics() {
    for (auto& cut : fresh_cuts_) {
        // Calc metrics
    }
    cuts_.insert(cuts_.end(), fresh_cuts_.begin(), fresh_cuts_.end());
}

Index Cutter::Select() {
    Index n_selected = 0;
    for (auto& cut : cuts_) {
        cut.selected = true;
        n_selected++;
    }
    return n_selected;
}

Index Cutter::Add() {
    Index m = model_.GetNCons();
    Index n = model_.GetNVars();

    Index n_added = 0;
    for (auto& cut : cuts_) {
        if (cut.selected) {
            // if (n_added >= 10) break;  // 10-й кат в бимва-п ломает
            // std::cout << "Adding " << cut;
            model_.PrepareConstraint(cut.lhs, {cut.rhs, kInf});
            n_added++;
        }
    }
    model_.Resize(m + n_added, n);
    model_.FinalizeAc();  // Todo можно без Srm2Scm, т.к. мы просто добавляем новые строки

    return n_added;
}

}  // namespace reshala
