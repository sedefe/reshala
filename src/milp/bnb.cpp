#include "reshala/milp/bnb.h"

namespace reshala {

Solution BnbSolver::Solve(const Solution& sol) {
    Node root(sol, model_.GetAllBounds());

    nodes.push_back(root);

    while (!nodes.empty()) {
        auto node = nodes.back();
        nodes.pop_back();

        if (mip_state_.Converged()) {
            break;
        }

        Index cand = branching.Branch(node.sol);

        Scalar x_floor = floor(node.sol.x[cand]);
        const Bounds& bnd = model_.GetBounds(cand);
        std::array<Bounds, 2> cand_bounds{{{bnd.le, x_floor}, {x_floor + 1, bnd.ri}}};

        model_.SetAllBounds(node.bounds);
        for (Index i = 0; i < 2; i++) {
            model_.SetBounds(cand, cand_bounds[i]);
            DualSimplex ds(model_);
            auto sol = ds.Solve();

            if (sol.status != LpStatus::kOptimal) {
                continue;
            }

            mip_state_.UpdDual(sol.y);
            if (model_.IsIntegerFeasible(sol.x)) {
                if (sol.y < best_sol_.y) {
                    best_sol_ = sol;
                }
                mip_state_.UpdPrimal(sol.y);
            } else {
                nodes.push_back(Node(sol, model_.GetAllBounds()));
            }
        }
    }

    return best_sol_;
}

}  // namespace reshala
