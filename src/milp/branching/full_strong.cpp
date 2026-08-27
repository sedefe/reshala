#include <array>
#include <limits>

#include "reshala/lp/dual_simplex.h"
#include "reshala/milp/branching.h"
#include "reshala/utils.h"

namespace reshala {

Index FullStrong::Branch(Node& parent, DualSimplex& ds) {
    Index candidate;
    bool candidate_used_ps;
    Scalar best_score;
    Index n_implied_bounds;
    Scalar gains[2];
    Solution sols[2];
    DsState ds_states[2];

    while (true) {
        candidate = -1;
        best_score = -kInf;
        n_implied_bounds = 0;

        // Родитель мог ухудшиться из-за подтягивания границ детьми
        if (parent.sol.y >= mip_tracker_.GetCutoff()) {
            parent.sol.status = LpStatus::kDropped;
            break;
        }

        for (Index iv = 0; iv < model_.GetNVars(); ++iv) {
            if (!model_.GetIntegrality(iv)) continue;
            const Scalar x_val = parent.sol.x[iv];
            if (IsZero(MinFraction(x_val))) continue;

            const Scalar floor_x = Floor(x_val);
            const Scalar ceil_x = floor_x + 1;

            Bounds orig_bnd = model_.GetBounds(iv);
            Bounds cand_bounds[2] = {{orig_bnd.le, floor_x}, {ceil_x, orig_bnd.ri}};
            Scalar dxs[2] = {x_val - floor_x, ceil_x - x_val};

            bool use_ps = false and hist_.IsEnough(iv);
            if (use_ps) {
                gains[0] = hist_.Estimate(iv, Direction::kLeft, dxs[0]);
                gains[1] = hist_.Estimate(iv, Direction::kRight, dxs[1]);
                sols[0].status = sols[1].status = LpStatus::kOptimal;
            } else {  // честно бранчим
                for (Index i = 0; i < 2; ++i) {
                    ds.Restore(parent.ds_state);
                    ds.SetBounds(iv, cand_bounds[i]);
                    sols[i] = ds.Solve(true, mip_tracker_.GetCutoff());
                    gains[i] = sols[i].y - parent.sol.y;

                    if (sols[i].status == LpStatus::kOptimal) {
                        hist_.Add(iv, Index2Dir(i), gains[i], dxs[i]);

                        // Катоф не прошёл => нахрен пошёл
                        if (sols[i].y >= mip_tracker_.GetCutoff()) {
                            sols[i].status = LpStatus::kDropped;
                            continue;
                        }

                        if (mip_tracker_.TestPrimal(sols[i])) {
                            ReportNewPrimal("FSB", sols[i].y);
                            if (mip_tracker_.Converged()) return 0;
                            // Хорошая была нода, но обрабатывать её дальше незачем
                            sols[i].status = LpStatus::kDropped;
                        } else {
                            ds_states[i] = ds.Store();
                            heur_manager_.Run(HeuristicTrigger::kFsb, model_, sols[i]);
                        }
                    }
                }

                // Если у кандидата нет детей, дропаем всю ноду
                if (sols[0].status != LpStatus::kOptimal and sols[1].status != LpStatus::kOptimal) {
                    parent.sol = InfeasibleSolution();
                    return 0;
                }

                // Если ровно один, второй занимает место родителя
                if (sols[0].status != LpStatus::kOptimal) {
                    n_implied_bounds++;
                    parent.domain.SetBounds(iv, cand_bounds[1]);
                    parent.ds_state = ds_states[1];
                    parent.sol = sols[1];
                    ds.SetBounds(iv, cand_bounds[1]);
                    continue;
                }
                if (sols[1].status != LpStatus::kOptimal) {
                    n_implied_bounds++;
                    parent.domain.SetBounds(iv, cand_bounds[0]);
                    parent.ds_state = ds_states[0];
                    parent.sol = sols[0];
                    ds.SetBounds(iv, cand_bounds[0]);
                    continue;
                }

                // Два => откатываем баунд
                ds.SetBounds(iv, orig_bnd);
            }

            Scalar score = (1.0 - kFsbMu) * std::min(gains[0], gains[1]) +
                           kFsbMu * std::max(gains[0], gains[1]);

            if (score > best_score) {
                best_score = score;
                candidate = iv;
                candidate_used_ps = use_ps;
            }
        }

        if (n_implied_bounds == 0) {
            break;
        }
    }

    if (candidate == -1) {  // Всё округлилось
        return 0;
    }

    // Это надо делать только тогда, когда два ребёнка
    const Scalar x_cand = parent.sol.x[candidate];
    const Scalar floor_cand = Floor(x_cand);
    const Scalar frac_cand = x_cand - floor_cand;
    const Bounds orig_bnd = model_.GetBounds(candidate);
    std::array<Bounds, 2> final_bounds{{{orig_bnd.le, floor_cand}, {floor_cand + 1, orig_bnd.ri}}};

    Index num_ch = 0;
    for (Index i = 0; i < 2; ++i) {
        ds.Restore(parent.ds_state);
        ds.SetBounds(candidate, final_bounds[i]);
        auto sol = ds.Solve(true, mip_tracker_.GetCutoff());
        heur_manager_.Run(HeuristicTrigger::kNode, model_, sol);
        children_[i] = Node(parent.level + 1, sol, model_.GetDomain(), ds.Store());

        if (candidate_used_ps) {
            if (sol.status == LpStatus::kOptimal) {
                hist_.Add(candidate, Index2Dir(i), sol.y - parent.sol.y,
                          i == 0 ? frac_cand : 1 - frac_cand);

                if (mip_tracker_.TestPrimal(sol)) {
                    ReportNewPrimal("Node", sol.y);
                }
            }
        }

        num_ch += (sol.status == LpStatus::kOptimal);
    }

    return num_ch;
}

}  // namespace reshala
