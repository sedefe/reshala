#include "reshala/milp/bnb.h"

namespace reshala {

std::ostream& operator<<(std::ostream& os, const BnbStats& stats) {
    os << "B&B: " << stats.n_nodes << " nodes visited, " << stats.n_dropped << " dropped\n";
    return os;
}

BnbSolver::BnbSolver(MilpModel& model, DualSimplex& ds, MipState& mip_state)
    : model_(model), ds_(ds), mip_state_(mip_state) {
    // branching_ = std::make_unique<FullStrong>(model);
    branching_ = std::make_unique<MostInfeasible>(model);
}

void BnbSolver::SolveRoot(Node& root) {
    FullStrongDomProp branching(model_);
    branching.Branch(root, ds_);
}

void BnbSolver::Solve(const Solution& relaxed) {
    Node root(1, relaxed, model_.GetDomain(), ds_.Store());
    SolveRoot(root);
    nodes.push_back(root);

    while (!nodes.empty()) {
        curr_node = std::move(nodes.back());
        nodes.pop_back();
        stats.n_nodes++;
        if (curr_node.sol.y >= mip_state_.GetCutoff()) {
            stats.n_dropped++;
            continue;
        }

        if (mip_state_.Converged()) {
            break;
        }

        model_.SetDomain(curr_node.domain);
        auto num_ch = branching_->Branch(curr_node, ds_);
        if (num_ch == 0) continue;

        Index best = branching_->FindBestChild();
        auto children_ids =
            (num_ch == 1) ? std::vector<Index>{best} : std::vector<Index>{1 - best, best};
        for (auto i : children_ids) {
            const Node& child = branching_->GetChild(i);
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

        UpdDual();
        DebugPrint();
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
    if (stats.n_nodes % 50 == 1) {
        std::cout
            << "===========================================================================\n";
        std::cout
            << "lev | LPiter | left       | right      | dual       | primal     | gap,%   \n";
        std::cout
            << "===========================================================================\n";
    }
    std::cout << FMT(3, 5) << curr_node.level << " | " << FMT(6, 5) << ds_.GetStats().n_iter
              << " | " << FMT(10, 5) << branching_->GetChild(0).sol.y << " | " << FMT(10, 5)
              << branching_->GetChild(1).sol.y << " | " << FMT(10, 5) << mip_state_.GetDual()
              << " | " << FMT(10, 5) << mip_state_.GetPrimal() << " | " << FMT(7, 4)
              << mip_state_.GetGap() * 1e2 << "\n"
              << FMT_DEFAULT;
}

}  // namespace reshala
