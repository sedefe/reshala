#include "reshala/milp/bnb.h"

namespace reshala {

std::ostream& operator<<(std::ostream& os, const BnbStats& stats) {
    os << "B&B: " << stats.n_nodes << " nodes visited, " << stats.n_dropped << " dropped\n";
    return os;
}

BnbSolver::BnbSolver(MilpModel& model, DualSimplex& ds, MipState& mip_state)
    : model_(model), ds_(ds), mip_state_(mip_state) {
    root_branching_ = std::make_unique<FullStrong>(model, mip_state);
    node_branching_ = std::make_unique<FullStrong>(model, mip_state);
    // node_branching_ = std::make_unique<MostInfeasible>(model, mip_state);
}

void BnbSolver::Solve(const Solution& relaxed) {
    Node root(1, relaxed, model_.GetDomain(), ds_.Store());
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
        std::unique_ptr<AbstractBranching>& branching_ =
            (curr_node.level == 1) ? root_branching_ : node_branching_;
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
            if (child.sol.y < mip_state_.GetCutoff()) {
                nodes.push_back(child);
            }
        }

        UpdDual();
        DebugPrint(branching_);
    }

    model_.SetDomain(root.domain);
}

void BnbSolver::UpdDual() {
    // В стеке могут оказаться только плохие ноды, поэтому праймал тоже считаем оценкой
    Scalar min_dual = mip_state_.GetPrimal();
    for (const auto& node : nodes) {
        if (node.sol.y < min_dual) {
            min_dual = node.sol.y;
        }
    }
    mip_state_.UpdDual(min_dual);
}

void BnbSolver::DebugPrint(std::unique_ptr<AbstractBranching>& branching) {
    if (stats.n_nodes % 50 == 1) {
        std::cout << "============================================================================="
                     "======\n";
        std::cout << "lev | LPiter | left         | right        | dual         | primal       |  "
                     "gap,%  \n";
        std::cout << "============================================================================="
                     "======\n";
    }
    std::cout << FMT(3, 5) << curr_node.level << " | " << FMT(6, 5)
              << FormatInteger(ds_.GetStats().n_iter) << " | " << FMT(12, 5)
              << branching->GetChild(0).sol.y << " | " << FMT(12, 5) << branching->GetChild(1).sol.y
              << " | " << FMT(12, 5) << mip_state_.GetDual() << " | " << FMT(12, 5)
              << mip_state_.GetPrimal() << " | " << FMT(7, 4) << mip_state_.GetGap() * 1e2 << "\n"
              << FMT_DEFAULT;
}

}  // namespace reshala
