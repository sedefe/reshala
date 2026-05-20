#include "reshala/presolve/presolve.h"

namespace reshala {

Presolver::Presolver(MilpModel& model) : model_(model), info_(model) {
    rules_.push_back(std::make_unique<Rule31>());
    rules_.push_back(std::make_unique<Rule32>());
    rules_.push_back(std::make_unique<Rule35>());
    rules_.push_back(std::make_unique<Rule41>());
    rules_.push_back(std::make_unique<Rule44>());
}

LpStatus Presolver::Presolve() {
    bool changed = true;
    Index pass = 0;
    Index max_passes = 5;

    while (changed && pass < max_passes && !info_.ProvenInfeasible()) {
        std::cout << "Presolve pass " << pass << ": " << model_.GetNCons() << " x "
                  << model_.GetNVars() << "\n";
        changed = false;
        info_.CalcActivities();

        for (auto& rule : rules_) {
            auto res = rule->Apply(info_, transforms_);

            if (info_.ProvenInfeasible()) break;
            if (info_.GetNDeletedCons() > 0) info_.CompressCons();
            if (info_.GetNDeletedVars() > 0) info_.CompressVars();

            if (res == RuleResult::kReduced) {
                changed = true;
            }
        }

        pass++;
    }
    std::cout << "After presolve: " << model_.GetNCons() << "  x " << model_.GetNVars() << "\n";

    if (info_.ProvenInfeasible()) {
        std::cout << "Presolve proved infeasibility\n";
        return LpStatus::kInfeasible;
    }
    if (model_.GetNVars() == 0) {
        assert(model_.GetNCons() == 0);
        std::cout << "Presolve reduced to empty\n";
        return LpStatus::kOptimal;
    }

    return LpStatus::kUnknown;
}

Solution Presolver::Postsolve(const Solution& sol) {
    if (sol.status != LpStatus::kOptimal) {
        return {sol.status, kInf, {}};
    }

    Solution res = sol;
    res.x.assign(info_.GetOrigNVars(), kNan);
    for (Index iv = 0; iv < sol.x.size(); iv++) {
        res.x[info_.GetOrigVarIdx()[iv]] = sol.x[iv];
    }

    for (auto it = transforms_.rbegin(); it != transforms_.rend(); ++it) {
        (*it)->Undo(res);
    }
    return res;
}

}  // namespace reshala
