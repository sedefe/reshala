#include <array>
#include <limits>

#include "reshala/lp/dual_simplex.h"
#include "reshala/milp/branching.h"
#include "reshala/utils.h"

namespace reshala {

Index FullStrongDomProp::Branch(Node& parent, DualSimplex& ds) {
    Index candidate = -1;
    Scalar best_score;

    bool feasible = true;
    Domain domain = parent.domain;
    while (feasible) {
        best_score = -kInf;
        Index n_implied_bounds = 0;
        for (Index iv = 0; iv < model_.GetNVars(); ++iv) {
            if (!model_.GetIntegrality(iv)) continue;
            const Scalar x_val = parent.sol.x[iv];
            const Scalar floor_x = Floor(x_val);
            const Scalar ceil_x = floor_x + 1;

            if (GetFraction(x_val) < kEpsZero) continue;
            if (candidate < 0) candidate = iv;

            Bounds bnd = model_.GetBounds(iv);
            std::array<Bounds, 2> cand_bounds{{{bnd.le, floor_x}, {ceil_x, bnd.ri}}};

            Scalar gains[2];
            LpStatus statuses[2];
            for (Index i = 0; i < 2; ++i) {
                ds.Restore(parent.ds_state);
                model_.SetBounds(iv, cand_bounds[i]);
                auto sol = ds.Solve(true);

                gains[i] = sol.y - parent.sol.y;
                statuses[i] = sol.status;
            }
            model_.SetBounds(iv, bnd);

            // Если у кандидата нет детей, дропаем всю ноду.
            if (statuses[0] != LpStatus::kOptimal and statuses[1] != LpStatus::kOptimal) {
                parent.sol = InfeasibleSolution();
                return 0;
            }
            if (statuses[0] != LpStatus::kOptimal) {
                n_implied_bounds++;
                domain.SetBounds(iv, {ceil_x, bnd.ri});
            }
            if (statuses[1] != LpStatus::kOptimal) {
                n_implied_bounds++;
                domain.SetBounds(iv, {bnd.le, floor_x});
            }

            Scalar score = (1.0 - kFsbMu) * std::min(gains[0], gains[1]) +
                           kFsbMu * std::max(gains[0], gains[1]);

            if ((candidate == -1) or (score > best_score)) {
                best_score = score;
                candidate = iv;
            }
        }
        if (n_implied_bounds > 0) {
            parent.domain = domain;
            ds.Restore(parent.ds_state);
            model_.SetDomain(parent.domain);

            parent.sol = ds.Solve(true);

            parent.ds_state = ds.Store();
        } else {
            break;
        }
    }

    // Prepare children
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
