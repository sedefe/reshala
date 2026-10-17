#include "reshala/presolve/presolve.h"

namespace reshala {

inline bool RuleSuccess(RuleResult status) {
    return status == RuleResult::kReduced or status == RuleResult::kInfeasible;
}

Presolver::Presolver(MilpModel& model) : model_(model) {
    row_rules_.push_back(std::make_unique<Rule31>());
    row_rules_.push_back(std::make_unique<Rule32>());
    row_rules_.push_back(std::make_unique<Rule33>());
    row_rules_.push_back(std::make_unique<Rule35>());
    row_rules_.push_back(std::make_unique<Rule36>());
    row_rules_.push_back(std::make_unique<Rule46>());

    col_rules_.push_back(std::make_unique<Rule41>());
    col_rules_.push_back(std::make_unique<Rule44>());
    col_rules_.push_back(std::make_unique<Rule48>());

    hard_rules_.push_back(std::make_unique<Rule47>());
    hard_rules_.push_back(std::make_unique<Rule52>());
    hard_rules_.push_back(std::make_unique<Rule76>());
    hard_rules_.push_back(std::make_unique<Rule72>());
    hard_rules_.push_back(std::make_unique<Rule53>());
}

RuleResult Presolver::SimplePresolve(bool verbose) {
    RuleResult sp_result = RuleResult::kUnchanged;

    std::vector<Index> indices;

    tracker_->TouchAll();
    tracker_->CalcLocks();
    tracker_->CalcActivities();

    bool changed = true;
    while (changed) {
        changed = false;

        // printf("start row presolve\n");
        tracker_->TouchAll();
        indices = tracker_->GetChangedCons().GetValues();
        tracker_->ClearChangedCons();
        for (auto& rule : row_rules_) {
            tracker_->CalcLocks();
            tracker_->CalcActivities();
            PresolveStat rule_stat = tracker_->stat;
            // printf("\t%s\n", rule->GetName().c_str());

            RuleResult rule_result = RuleResult::kUnchanged;
            for (Index ic : indices) {
                if (tracker_->ConIsDeleted(ic)) continue;
                RuleResult row_result = rule->Apply(*tracker_, ic);

                if (RuleSuccess(row_result)) rule_result = row_result;
                if (row_result == RuleResult::kInfeasible) break;
            }
            // printf("\t\t%d\n", rule_status);

            if (RuleSuccess(rule_result)) {
                sp_result = rule_result;
                changed = true;

                if (rule_result == RuleResult::kReduced) {
                    rule_stat = tracker_->stat - rule_stat;
                    if (verbose) PrintStat(*rule, rule_stat);
                }
                if (rule_result == RuleResult::kInfeasible) {
                    break;
                }
            }
        }
        if (sp_result == RuleResult::kInfeasible) break;
        if (tracker_->GetNDeletedCons() > 0) tracker_->CompressCons();
        if (tracker_->GetNDeletedVars() > 0) tracker_->CompressVars();
        tracker_->CalcLocks();
        tracker_->CalcActivities();

        // printf("start col presolve\n");
        indices = tracker_->GetChangedVars().GetValues();
        tracker_->ClearChangedVars();
        for (auto& rule : col_rules_) {
            tracker_->CalcLocks();
            tracker_->CalcActivities();
            PresolveStat rule_stat = tracker_->stat;
            // printf("\t%s\n", rule->GetName().c_str());

            RuleResult rule_result = RuleResult::kUnchanged;
            for (Index iv : indices) {
                if (tracker_->VarIsDeleted(iv)) continue;
                RuleResult col_result = rule->Apply(*tracker_, iv);

                if (RuleSuccess(col_result)) rule_result = col_result;
                if (col_result == RuleResult::kInfeasible) break;
            }
            // printf("\t\t%d\n", rule_status);

            if (RuleSuccess(rule_result)) {
                sp_result = rule_result;
                changed = true;

                if (rule_result == RuleResult::kReduced) {
                    rule_stat = tracker_->stat - rule_stat;
                    if (verbose) PrintStat(*rule, rule_stat);
                }
                if (rule_result == RuleResult::kInfeasible) {
                    break;
                }
            }
        }
        if (sp_result == RuleResult::kInfeasible) break;
        if (tracker_->GetNDeletedCons() > 0) tracker_->CompressCons();
        if (tracker_->GetNDeletedVars() > 0) tracker_->CompressVars();
        tracker_->CalcLocks();
        tracker_->CalcActivities();
    }

    return sp_result;
}

LpStatus Presolver::Presolve(bool verbose, bool hard) {
    tracker_ = std::make_unique<ModelTracker>(model_);

    if (verbose) PrintHeader();
    status_ = SimplePresolve(verbose);
    // exit(0);

    if (hard and !Solved()) {
        bool changed = true;
        while (changed and !Solved()) {
            changed = false;
            for (auto& rule : hard_rules_) {
                PresolveStat rule_stat = tracker_->stat;
                // printf("\t%s\n", rule->GetName().c_str());
                RuleResult rule_result = rule->Apply(*tracker_);

                if (RuleSuccess(rule_result)) {
                    status_ = rule_result;
                    changed = true;

                    if (rule_result == RuleResult::kReduced) {
                        rule_stat = tracker_->stat - rule_stat;
                        if (verbose) PrintStat(*rule, rule_stat);
                        status_ = SimplePresolve(verbose);
                    }
                    break;
                }
            }
        }
    }

    // if (verbose) std::cout << "After presolve: " << model_.StatString() << "\n";
    if (verbose) printf("%d After presolve: %s\n", status_, model_.StatString().c_str());
    status_ = RuleResult::kInfeasible;

    if (status_ == RuleResult::kInfeasible) {
        if (verbose) std::cout << "Presolve proved infeasibility\n";
        return LpStatus::kInfeasible;
    }
    if (model_.GetNVars() == 0) {
        assert(model_.GetNCons() == 0);
        if (verbose) std::cout << "Presolve reduced to empty\n";
        return LpStatus::kOptimal;
    }

    {  // Check if there are free vars
        Index n_free = 0;
        for (Index iv = 0; iv < model_.GetNVars(); iv++)
            n_free += model_.GetType(iv) == BndType::kFree;
        if (n_free > 0) {
            std::cerr << "No support for free variables yet, but there are " << n_free
                      << " of them\n";
            return LpStatus::kError;
        }
    }

    return LpStatus::kUnknown;
}

Solution Presolver::Postsolve(const Solution& sol) {
    switch (sol.status) {
        case LpStatus::kInfeasible:
        case LpStatus::kDropped:
            return Solution::Infeasible();
            break;
        case LpStatus::kError:
            return {sol.status, kNan, {}};
            break;
        default:
            break;
    }

    assert(sol.status == LpStatus::kOptimal);

    Solution res = sol;
    res.y = model_.GetObj().evaluate(sol.x);

    res.x.assign(tracker_->GetOrigNVars(), kNan);
    for (Index iv = 0; iv < sol.x.size(); iv++) {
        res.x[tracker_->GetOrigVarIdx()[iv]] = sol.x[iv];
    }

    for (auto it = tracker_->GetTransforms().rbegin(); it != tracker_->GetTransforms().rend();
         ++it) {
        (*it)->Undo(res);
    }
    return res;
}

void Presolver::PrintHeader() const {
    std::cout << "Rule name       RmC   RmV   ChB   ChR   ChC    Nnz\n";
}

void Presolver::PrintStat(const Rule& rule, const PresolveStat& stat) const {
    std::cout << std::left << std::setw(12) << rule.GetName() << ": " << std::right << std::setw(5)
              << stat.n_rm_con << " " << std::setw(5) << stat.n_rm_var << " " << std::setw(5)
              << stat.n_ch_bnd << " " << std::setw(5) << stat.n_ch_rhs << " " << std::setw(5)
              << stat.n_ch_coeff << " " << std::setw(6) << model_.GetNnz() << "\n";
}

}  // namespace reshala
