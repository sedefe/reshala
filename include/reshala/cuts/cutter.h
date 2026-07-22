#pragma once

#include "reshala/cuts/config.h"
#include "reshala/cuts/cut.h"
#include "reshala/cuts/generators/cmir.h"
#include "reshala/cuts/generators/probing.h"
#include "reshala/milp/utils.h"

namespace reshala {

struct CgStats {
    Index n_generated = 0;
    Index n_active = 0;
};
struct CutterStats {
    CutterStats() {
        cg_stats.emplace(CutType::kProbing, CgStats());
        cg_stats.emplace(CutType::kCmir, CgStats());
    }
    std::unordered_map<CutType, CgStats> cg_stats;
};
std::ostream& operator<<(std::ostream& os, const CutterStats& stats);

class Cutter {
   public:
    Cutter(MilpModel& model, const Presolver& presolver, DualSimplex& ds, MipState& mip_state);

    void Run(Solution& sol);

    inline const CutterStats& GetStats() const { return stats; }

   private:
    MilpModel& model_;
    const Presolver& presolver_;
    DualSimplex& ds_;
    MipState& mip_state_;
    Solution sol_;

    CutterStats stats;

    Index n_round_ = 0;

    Index max_cuts_;
    Index max_support_;

    std::vector<Cut> pool_;

    Index Generate(const Solution& sol);
    void Estimate();
    void Filter();
    Index Select();
    Index Add();

    bool AnyNewSelected();
};

}  // namespace reshala
