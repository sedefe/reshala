#include "reshala/presolve/presolve.h"

namespace reshala {

inline bool RuleSuccess(RuleResult status) {
    return status == RuleResult::kReduced or status == RuleResult::kInfeasible;
}

Presolver::Presolver(MilpModel& model) : model_(model) {
    simple_rules_.push_back(std::make_unique<Rule31>());
    simple_rules_.push_back(std::make_unique<Rule32>());
    simple_rules_.push_back(std::make_unique<Rule33>());
    simple_rules_.push_back(std::make_unique<Rule35>());
    simple_rules_.push_back(std::make_unique<Rule36>());
    simple_rules_.push_back(std::make_unique<Rule41>());
    simple_rules_.push_back(std::make_unique<Rule44>());
    simple_rules_.push_back(std::make_unique<Rule46>());
    simple_rules_.push_back(std::make_unique<Rule48>());
    hard_rules_.push_back(std::make_unique<Rule47>());
    hard_rules_.push_back(std::make_unique<Rule52>());
    hard_rules_.push_back(std::make_unique<Rule76>());
    hard_rules_.push_back(std::make_unique<Rule72>());
    hard_rules_.push_back(std::make_unique<Rule53>());
}

void Presolver::SimplePresolve(bool verbose) {
    status_ = RuleResult::kUnchanged;
    RuleResult rule_result;

    bool changed = true;
    while (changed and status_ != RuleResult::kInfeasible) {
        changed = false;

        for (auto& rule : simple_rules_) {
            PresolveStat rule_stat = tracker_->stat;
            rule_result = rule->Apply(*tracker_);

            if (RuleSuccess(rule_result)) {
                status_ = rule_result;
                changed = true;
                rule_stat = tracker_->stat - rule_stat;
                if (verbose) PrintStat(*rule, rule_stat);
            }
            if (rule_result == RuleResult::kInfeasible) break;

            if (tracker_->GetNDeletedCons() > 0) tracker_->CompressCons();
            if (tracker_->GetNDeletedVars() > 0) tracker_->CompressVars();

            {
                auto locks = tracker_->GetLocks();
                tracker_->CalcLocks();
                int nz = 0;
                if (verbose) {
                    for (Index i = 0; i < model_.GetNVars(); i++) {
                        auto ld1 = locks.GetNLocks(i, LockType::kDown);
                        auto lu1 = locks.GetNLocks(i, LockType::kUp);
                        auto lu2 = tracker_->GetLocks().GetNLocks(i, LockType::kUp);
                        auto ld2 = tracker_->GetLocks().GetNLocks(i, LockType::kDown);
                        if (lu1 != lu2 or ld1 != ld2) {
                            nz++;
                            printf("Zalupa x%-4d: [%2d %2d] vs [%2d %2d]\n", i, ld1, lu1, ld2, lu2);
                        }
                    }
                }
                if (verbose and nz > 0) {
                    printf("%d zalup\n", nz);
                    exit(0);
                }
            }
        }
    }
}

LpStatus Presolver::Presolve(bool verbose, bool hard) {
    tracker_ = std::make_unique<ModelTracker>(model_);

    if (verbose) PrintHeader();

    SimplePresolve(verbose);

    if (hard and status_ != RuleResult::kInfeasible) {
        bool changed = true;
        while (changed and status_ != RuleResult::kInfeasible) {
            changed = false;
            for (auto& rule : hard_rules_) {
                PresolveStat rule_stat = tracker_->stat;
                RuleResult rule_result = rule->Apply(*tracker_);

                if (RuleSuccess(rule_result)) {
                    status_ = rule_result;
                    changed = true;
                    rule_stat = tracker_->stat - rule_stat;
                    if (verbose) PrintStat(*rule, rule_stat);

                    if (rule_result == RuleResult::kReduced) {
                        if (tracker_->GetNDeletedCons() > 0) tracker_->CompressCons();
                        if (tracker_->GetNDeletedVars() > 0) tracker_->CompressVars();
                        SimplePresolve(verbose);
                    }
                    break;
                }
            }
        }
    }

    if (verbose) std::cout << "After presolve: " << model_.StatString() << "\n";

    if (status_ == RuleResult::kInfeasible) {
        if (verbose) std::cout << "Presolve proved infeasibility\n";
        return LpStatus::kInfeasible;
    }
    if (model_.GetNVars() == 0) {
        assert(model_.GetNCons() == 0);
        if (verbose) std::cout << "Presolve reduced to empty\n";
        return LpStatus::kOptimal;
    }

    auto locks = tracker_->GetLocks();
    tracker_->CalcLocks();
    int nz = 0;
    if (verbose) {
        for (Index i = 0; i < model_.GetNVars(); i++) {
            auto ld1 = locks.GetNLocks(i, LockType::kDown);
            auto lu1 = locks.GetNLocks(i, LockType::kUp);
            auto lu2 = tracker_->GetLocks().GetNLocks(i, LockType::kUp);
            auto ld2 = tracker_->GetLocks().GetNLocks(i, LockType::kDown);
            if (lu1 != lu2 or ld1 != ld2) {
                nz++;
                printf("Zalupa x%-4d: [%2d %2d] vs [%2d %2d]\n", i, ld1, lu1, ld2, lu2);
            }
        }
    }
    if (verbose) printf("%d zalup\n", nz);

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
