#include "reshala/milp/bnb.h"

namespace reshala {

Solution BnbSolver::Solve(const Solution& sol) {
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
                    printf("New integer solution: %8.5g\n", child.sol.y);
                }
            } else {
                if (child.sol.y < mip_state_.GetCutoff()) {
                    nodes.push_back(child);
                }
            }
        }
    }

    return mip_state_.GetBestSol();
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
        printf("=========================================================\n");
        printf("lev | left     | right    | dual     | primal   | gap    \n");
        printf("=========================================================\n");
    }
    printf("%3ld | %8.5g | %8.5g | %8.5g | %8.5g | %6.2f%%\n", nodes.size(),
           branching.GetChild(0).sol.y, branching.GetChild(1).sol.y, mip_state_.GetDual(),
           mip_state_.GetPrimal(), mip_state_.GetGap() * 1e2);
}

}  // namespace reshala
