#include <array>
#include <limits>

#include "reshala/lp/dual_simplex.h"
#include "reshala/milp/branching.h"
#include "reshala/utils.h"

namespace reshala {

Index FullStrong::Branch(Node& parent, DualSimplex& ds) {
    Index candidate;
    Scalar best_score;
    Index n_implied_bounds;

    while (true) {
        candidate = -1;
        best_score = -kInf;
        n_implied_bounds = 0;

        for (Index iv = 0; iv < model_.GetNVars(); ++iv) {
            if (!model_.GetIntegrality(iv)) continue;
            const Scalar x_val = parent.sol.x[iv];
            if (GetFraction(x_val) <= kEpsZero) continue;

            const Scalar floor_x = Floor(x_val);
            const Scalar ceil_x = floor_x + 1;

            Bounds orig_bnd = model_.GetBounds(iv);
            std::array<Bounds, 2> cand_bounds{{{orig_bnd.le, floor_x}, {ceil_x, orig_bnd.ri}}};

            Solution sols[2];
            DsState ds_states[2];
            for (Index i = 0; i < 2; ++i) {
                ds.Restore(parent.ds_state);
                model_.SetBounds(iv, cand_bounds[i]);
                sols[i] = ds.Solve(true);

                // Катоф не прошёл => нахрен пошёл
                if (sols[i].y > mip_state_.GetCutoff()) {
                    sols[i].status = LpStatus::kInfeasible;
                }

                if (sols[i].status == LpStatus::kOptimal) {
                    if (mip_state_.TestPrimal(sols[i])) {
                        std::cout << "FSB: New integer solution: " << FMT(10, 5) << sols[i].y
                                  << "\n";
                        if (mip_state_.Converged()) return 0;
                    }
                    ds_states[i] = ds.Store();
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
                model_.SetBounds(iv, cand_bounds[1]);
                continue;
            }
            if (sols[1].status != LpStatus::kOptimal) {
                n_implied_bounds++;
                parent.domain.SetBounds(iv, cand_bounds[0]);
                parent.ds_state = ds_states[0];
                parent.sol = sols[0];
                model_.SetBounds(iv, cand_bounds[0]);
                continue;
            }

            // Два => откатываем баунд и считаем скор
            model_.SetBounds(iv, orig_bnd);

            Scalar gains[2] = {sols[0].y - parent.sol.y, sols[1].y - parent.sol.y};
            Scalar score = (1.0 - kFsbMu) * std::min(gains[0], gains[1]) +
                           kFsbMu * std::max(gains[0], gains[1]);

            if (score > best_score) {
                best_score = score;
                candidate = iv;
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
        model_.SetBounds(candidate, final_bounds[i]);
        auto sol = ds.Solve(true);
        children_[i] = Node(parent.level + 1, sol, model_.GetDomain(), ds.Store());

        num_ch += (sol.status == LpStatus::kOptimal);
    }

    return num_ch;
}

}  // namespace reshala
