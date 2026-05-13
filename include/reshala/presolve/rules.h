#pragma once

#include <memory>
#include <string>

#include "reshala/presolve/transforms.h"

namespace reshala {

enum class RuleType { kTrivial, kFast, kMedium, kExhaustive, kUnknown };
enum class RuleResult { kSkipped, kUnchanged, kReduced, kInfeasible, kUnknown };

class Rule {
   public:
    virtual ~Rule() = default;

    virtual RuleResult Apply(MilpModel& model,
                             std::vector<std::unique_ptr<Transform>>& transforms_) = 0;

    virtual std::string Name() const = 0;
};

class Rule41 : public Rule {
    RuleResult Apply(MilpModel& model, std::vector<std::unique_ptr<Transform>>& transforms_);
    std::string Name() const { return "4.1 Removal of fixed variables"; }
};

}  // namespace reshala
