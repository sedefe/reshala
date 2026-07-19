#include "reshala/cuts/cutter.h"

#include "reshala/lina/core/operators.h"

namespace reshala {

bool CutCompare(const Cut& c1, const Cut& c2) {
    if (c1.removed != c2.removed) {
        return c1.removed < c2.removed;  // удалённый кат всегда хуже
    }
    return c1.quality > c2.quality;
}

Cutter::Cutter(MilpModel& model, const Presolver& presolver, DualSimplex& ds, MipState& mip_state)
    : model_(model), presolver_(presolver), ds_(ds), mip_state_(mip_state) {
    auto m = model.GetNCons();
    auto n = model.GetNVars();

    max_cuts_ = Index(kMaxCutsFactor * m);
    max_support_ = std::max(2, Index(kMaxRelSupport * n));  // To guarantee probing is allowed
}

void Cutter::Run(Solution& sol) {
    for (n_round_ = 0; n_round_ < kMaxRounds; n_round_++) {
        sol_ = sol;

        auto n_generated = Generate(sol);
        Estimate();
        Filter();
        auto n_selected = Select();
        if (!AnyNewSelected()) break;

        auto n_added = Add();
        std::cout << "Added " << n_added << " cuts, size: " << model_.GetNCons() << " x "
                  << model_.GetNVars() << "\n";

        if (n_added > 0) {
            ds_.SetModel(model_);
            sol = ds_.Solve(false);
            if (sol.y > mip_state_.GetDual()) {
                mip_state_.UpdDual(sol.y);
            }
        }

        if (sol.status != LpStatus::kOptimal) break;
    }

    std::cout << "After cutter: " << sol.y << "\n";
}

Index Cutter::Generate(const Solution& sol) {
    Index prev_size = pool_.size();

    std::vector<std::unique_ptr<AbstractCg>> generators_;
    generators_.push_back(std::make_unique<ProbingCg>(model_, presolver_, ds_));

    for (auto& cg : generators_) {
        Index k = pool_.size();
        cg->Generate(sol, pool_);
        k = pool_.size() - k;
        std::cout << cg->GetName() << " generated " << k << " cuts\n";
    }

    for (auto i = prev_size; i < pool_.size(); i++) {
        pool_[i].round = n_round_;
    }

    return pool_.size() - prev_size;
}

void Cutter::Estimate() {
    for (auto& cut : pool_) {
        if (!cut.removed) {
            cut.CalcMetrics(sol_, model_);
            cut.quality = cut.m_obj_parallelism;
        }
    }
    std::sort(pool_.begin(), pool_.end(), CutCompare);
}

void Cutter::Filter() {
    {  // Filter by support
        Index n_filtered = 0;
        for (auto& cut : pool_) {
            if (cut.removed) continue;
            if (cut.m_support > max_support_) {
                cut.removed = true;
                n_filtered++;
            }
        }
        // std::cout << "Filtered " << n_filtered << " by support\n";
    }

    {  // Filter by age
        Index n_filtered = 0;
        for (auto& cut : pool_) {
            if (cut.removed) continue;
            cut.age++;
            if (cut.age > kMaxAge) {
                cut.removed = true;
                n_filtered++;
            }
        }
        // std::cout << "Filtered " << n_filtered << " by age\n";
    }

    {  // Filter parallel cuts
        Index n_filtered = 0;
        for (Index ic1 = 0; ic1 < Index(pool_.size()); ic1++) {
            if (pool_[ic1].removed) continue;
            for (Index ic2 = ic1 + 1; ic2 < Index(pool_.size()); ic2++) {
                if (pool_[ic2].removed) continue;
                auto cosine2 = Cos2(pool_[ic1].lhs, pool_[ic2].lhs);
                if (cosine2 > kThdCos2) {
                    pool_[ic2].removed = true;
                    n_filtered++;
                }
            }
        }
        // std::cout << "Filtered " << n_filtered << " parallel\n";
    }
}

Index Cutter::Select() {
    Index n_selected = 0;
    auto best_quality = -kInf;
    for (const auto& cut : pool_) {
        if (!cut.removed) {
            best_quality = cut.quality;
            break;
        }
    }
    for (auto& cut : pool_) {
        cut.selected = ((!cut.removed) & (cut.quality >= kMinRelQuality * best_quality));
        if (cut.selected) {
            cut.age = 1;
            n_selected++;
        }
        if (n_selected >= max_cuts_) break;
    }

    return n_selected;
}

Index Cutter::Add() {
    Index m = model_.GetNCons();
    Index n = model_.GetNVars();

    Index n_added = 0;
    for (auto& cut : pool_) {
        if (cut.selected) {
            model_.PrepareConstraint(cut.lhs, {cut.rhs, kInf});
            n_added++;
        }
    }
    model_.Resize(m + n_added, n);
    model_.FinalizeAc();  // Todo можно без Srm2Scm, т.к. мы просто добавляем новые строки

    return n_added;
}

bool Cutter::AnyNewSelected() {
    for (const auto& cut : pool_) {
        if (cut.round == n_round_ and cut.selected) return true;
    }
    return false;
}

}  // namespace reshala
