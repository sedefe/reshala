#include "reshala/presolve/presolve.h"

namespace reshala {

Presolver::Presolver(MilpModel& model) : model_(model) {
    rules_.push_back(std::make_unique<Rule41>());
}

void Presolver::Presolve() {
    bool changed = true;
    Index pass = 0;
    Index max_passes = 5;

    while (changed && pass < max_passes) {
        changed = false;
        for (auto& rule : rules_) {
            if (rule->Apply(model_, transforms_) == RuleResult::kReduced) {
                changed = true;
            }
        }
        pass++;
    }
}

Solution Presolver::Postsolve(const Solution& sol) {
    Solution res = sol;
    for (auto it = transforms_.rbegin(); it != transforms_.rend(); ++it) {
        (*it)->Undo(res);
    }
    return res;
}

}  // namespace reshala
