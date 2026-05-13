#include "reshala/milp/bnb.h"

namespace reshala {

void BnbSolver::Solve(const Solution& sol) {
    Node root(1, sol, model_.GetDomain());
    nodes.push_back(root);

    while (!nodes.empty()) {
        UpdDual();

        Node node = nodes.back();
        nodes.pop_back();
        n_nodes_++;
        model_.SetDomain(node.domain);

        if (mip_state_.Converged()) {
            break;
        }

        Index candidate = branching.Branch(node);
        DebugPrint();

        Index best = branching.FindBestChild();
        auto children_ids = {1 - best, best};
        for (auto i : children_ids) {
            const Node& child = branching.GetChild(i);
            if (child.sol.status != LpStatus::kOptimal) {
                // Todo: run conflict analysis
                continue;
            }
            if (model_.IsIntegerFeasible(child.sol.x)) {
                if (mip_state_.TestPrimal(child.sol)) {
                    std::cout << "New integer solution: " << child.sol.y << "\n";
                }
            } else {
                if (child.sol.y < mip_state_.GetCutoff()) {
                    nodes.push_back(child);
                }
            }
        }
    }

    model_.SetDomain(root.domain);
}

void BnbSolver::UpdDual() {
    Scalar min_dual = kInf;
    for (const auto& node : nodes) {
        if (node.sol.y < min_dual) {
            min_dual = node.sol.y;
        }
    }
    mip_state_.UpdDual(min_dual);
}

void BnbSolver::DebugPrint() {
    if (n_nodes_ % 50 == 1) {
        std::cout << "==================================================================\n";
        std::cout << "lev | left       | right      | dual       | primal     | gap,%   \n";
        std::cout << "==================================================================\n";
    }
    std::cout << FMT(3, 5) << nodes.size() << " | " << FMT(10, 5) << branching.GetChild(0).sol.y
              << " | " << FMT(10, 5) << branching.GetChild(1).sol.y << " | " << FMT(10, 5)
              << mip_state_.GetDual() << " | " << FMT(10, 5) << mip_state_.GetPrimal() << " | "
              << FMT(7, 4) << mip_state_.GetGap() * 1e2 << "\n"
              << FMT_DEFAULT;
}

}  // namespace reshala
