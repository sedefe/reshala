#include "reshala/presolve/presolve.h"

namespace reshala {

Presolver::Presolver(MilpModel& model) : model_(model), tracker_(model) {
    std::vector<std::unique_ptr<Rule>> rules_;
    rules_.push_back(std::make_unique<Rule31>(RuleType::kFast));
    rules_.push_back(std::make_unique<Rule32>(RuleType::kMedium));
    rules_.push_back(std::make_unique<Rule33>(RuleType::kFast));
    rules_.push_back(std::make_unique<Rule35>(RuleType::kFast));
    rules_.push_back(std::make_unique<Rule36>(RuleType::kFast));
    rules_.push_back(std::make_unique<Rule41>(RuleType::kFast));
    rules_.push_back(std::make_unique<Rule44>(RuleType::kFast));
    rules_.push_back(std::make_unique<Rule52>(RuleType::kMedium));

    for (auto& rule : rules_) {
        RuleType type = rule->type;
        rule_map_[type].push_back(std::move(rule));
    }
}

LpStatus Presolver::Presolve() {
    RuleType curr_level = RuleType::kFast;
    Index pass = 0;

    tracker_.CalcActivities();

    PrintHeader();
    while (!tracker_.ProvenInfeasible()) {
        bool changed = false;
        PresolveStat round_stat = tracker_.stat;
        for (auto& rule : rule_map_[curr_level]) {
            PresolveStat rule_stat = tracker_.stat;
            auto res = rule->Apply(tracker_);
            if (tracker_.ProvenInfeasible()) break;
            if (res == RuleResult::kReduced) {
                changed = true;
                rule_stat = tracker_.stat - rule_stat;
                PrintStat(*rule, rule_stat);
            }
        }
        if (tracker_.GetNDeletedCons() > 0) tracker_.CompressCons();
        if (tracker_.GetNDeletedVars() > 0) tracker_.CompressVars();

        round_stat = tracker_.stat - round_stat;
        if (!changed) {
            curr_level = NextLevel(curr_level);
            if (curr_level == RuleType::kUnknown) break;
        } else {
            curr_level = RuleType::kFast;
        }
    }

    std::cout << "After presolve: " << model_.GetNCons() << " x " << model_.GetNVars() << ", "
              << model_.GetNnz() << " nnz\n";

    if (tracker_.ProvenInfeasible()) {
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
    res.x.assign(tracker_.GetOrigNVars(), kNan);
    for (Index iv = 0; iv < sol.x.size(); iv++) {
        res.x[tracker_.GetOrigVarIdx()[iv]] = sol.x[iv];
    }

    for (auto it = tracker_.GetTransforms().rbegin(); it != tracker_.GetTransforms().rend(); ++it) {
        (*it)->Undo(res);
    }
    return res;
}

void Presolver::PrintHeader() const {
    std::cout << "Rule name       RmC   RmV   ChB   ChR   ChC    Nnz\n";
}

void Presolver::PrintStat(const Rule& rule, const PresolveStat& stat) const {
    std::cout << std::left << std::setw(12) << rule.Name() << ": " << std::right << std::setw(5)
              << stat.n_rm_con << " " << std::setw(5) << stat.n_rm_var << " " << std::setw(5)
              << stat.n_ch_bnd << " " << std::setw(5) << stat.n_ch_rhs << " " << std::setw(5)
              << stat.n_ch_coeff << " " << std::setw(6) << model_.GetNnz() << "\n";
}

}  // namespace reshala
