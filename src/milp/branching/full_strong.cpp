#include <array>
#include <limits>

#include "reshala/lp/dual_simplex.h"
#include "reshala/milp/branching.h"
#include "reshala/utils.h"

namespace reshala {

Index FullStrong::Branch(const Node& parent, DualSimplex& ds) {
    const Scalar mu = 0.5;

    Index candidate = -1;
    bool have_1ch_cand = false;
    Scalar best_score = -kInf;

    for (Index iv = 0; iv < model_.GetNVars(); ++iv) {
        if (!model_.GetIntegrality(iv)) continue;

        const Scalar x_val = parent.sol.x[iv];
        const Scalar floor_x = Floor(x_val);
        const Scalar ceil_x = floor_x + 1;

        if (GetFraction(x_val) < kEpsZero) continue;

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
            candidate = iv;
            break;
        }

        Scalar score;
        // Кандидаты с одним ребёнком в приоритете.
        // Среди таких выбран будет тот, у кого выше обжектив единственного ребёнка.
        if (statuses[0] != LpStatus::kOptimal or statuses[1] != LpStatus::kOptimal) {
            if (statuses[0] == LpStatus::kOptimal) score = gains[0];
            if (statuses[1] == LpStatus::kOptimal) score = gains[1];

            if ((!have_1ch_cand) or (score > best_score)) {
                best_score = score;
                candidate = iv;
            }
            have_1ch_cand = true;
        }
        if (have_1ch_cand) continue;

        // Остаются кандидаты с двумя детьми.
        score =
            (1.0 - kFsbMu) * std::min(gains[0], gains[1]) + kFsbMu * std::max(gains[0], gains[1]);

        if ((candidate == -1) or (score > best_score)) {
            best_score = score;
            candidate = iv;
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
    model_.SetBounds(candidate, orig_bnd);

    return num_ch;
}

}  // namespace reshala
