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

    while (true) {
        candidate = -1;
        best_score = -kInf;
        n_implied_bounds = 0;

        for (Index iv = 0; iv < model_.GetNVars(); ++iv) {
            if (!model_.GetIntegrality(iv)) continue;
            const Scalar x_val = parent.sol.x[iv];
            if (MinFraction(x_val) <= kEpsZero) continue;

            const Scalar floor_x = Floor(x_val);
            const Scalar ceil_x = floor_x + 1;

            Bounds orig_bnd = model_.GetBounds(iv);
            Bounds cand_bounds[2] = {{orig_bnd.le, floor_x}, {ceil_x, orig_bnd.ri}};
            Scalar dxs[2] = {x_val - floor_x, ceil_x - x_val};

            bool use_ps = false and hist_.IsEnough(iv);

            Solution sols[2];
            DsState ds_states[2];
            for (Index i = 0; i < 2; ++i) {
                Direction dir = Index2Dir(i);
                if (use_ps) {
                    gains[i] = hist_.Estimate(iv, dir, dxs[i]);
                    sols[i].status = LpStatus::kOptimal;
                    continue;
                }

                ds.Restore(parent.ds_state);
                ds.SetBounds(iv, cand_bounds[i]);
                sols[i] = ds.Solve(true);
                gains[i] = sols[i].y - parent.sol.y;

                if (sols[i].status != LpStatus::kInfeasible) {
                    hist_.Add(iv, dir, gains[i], dxs[i]);
                }

                // Катоф не прошёл => нахрен пошёл
                if (sols[i].y >= mip_state_.GetCutoff()) {
                    sols[i].status = LpStatus::kDropped;
                    continue;
                }

                if (sols[i].status == LpStatus::kOptimal) {
                    if (mip_state_.TestPrimal(sols[i])) {
                        std::cout << "FSB: New integer solution: " << FMT(10, 5) << sols[i].y
                                  << "\n";
                        if (mip_state_.Converged()) return 0;
                        // Хорошая была нода, но обрабатывать её дальше незачем
                        sols[i].status = LpStatus::kDropped;
                    } else {
                        ds_states[i] = ds.Store();
                    }
                }
            }

            if (!use_ps) {
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
    const Bounds orig_bnd = model_.GetBounds(candidate);
    std::array<Bounds, 2> final_bounds{{{orig_bnd.le, floor_cand}, {floor_cand + 1, orig_bnd.ri}}};

    Index num_ch = 0;
    for (Index i = 0; i < 2; ++i) {
        ds.Restore(parent.ds_state);
        ds.SetBounds(candidate, final_bounds[i]);
        auto sol = ds.Solve(true);
        children_[i] = Node(parent.level + 1, sol, model_.GetDomain(), ds.Store());

        Scalar f = Fraction(parent.sol.x[candidate]);
        if (candidate_used_ps and sol.status != LpStatus::kInfeasible) {
            hist_.Add(candidate, Index2Dir(i), sol.y - parent.sol.y, i == 0 ? f : 1 - f);
        }

        num_ch += (sol.status == LpStatus::kOptimal);
    }

    return num_ch;
}

}  // namespace reshala
