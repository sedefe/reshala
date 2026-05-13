#pragma once

#include <memory>
#include <string>

#include "reshala/presolve/transforms.h"

namespace reshala {

enum class RuleType { kTrivial, kFast, kMedium, kExhaustive, kUnknown };
enum class RuleResult { kSkipped, kUnchanged, kReduced, kInfeasible, kUnknown };

class Rule {
   public:
    Rule(RuleType t) : type(t) {}
    virtual ~Rule() = default;

    virtual RuleResult Apply(ModelInfo& info,
                             std::vector<std::unique_ptr<Transform>>& transforms_) = 0;

    virtual std::string Name() const = 0;
    RuleType type;  // New int field
};

class Rule31 : public Rule {
   public:
    Rule31() : Rule(RuleType::kFast) {}
    RuleResult Apply(ModelInfo& info, std::vector<std::unique_ptr<Transform>>& transforms_);
    std::string Name() const { return "3.1 Removal of redundant constraints"; }
};

class Rule32 : public Rule {
   public:
    Rule32() : Rule(RuleType::kFast) {}
    RuleResult Apply(ModelInfo& info, std::vector<std::unique_ptr<Transform>>& transforms_);
    std::string Name() const { return "3.2 Bound strengthening"; }
};

class Rule41 : public Rule {
   public:
    Rule41() : Rule(RuleType::kFast) {}
    RuleResult Apply(ModelInfo& info, std::vector<std::unique_ptr<Transform>>& transforms_);
    std::string Name() const { return "4.1 Removal of fixed variables"; }
};

class Rule44 : public Rule {
   public:
    Rule44() : Rule(RuleType::kFast) {}
    RuleResult Apply(ModelInfo& info, std::vector<std::unique_ptr<Transform>>& transforms_);
    std::string Name() const { return "4.4 Dual fixing"; }
};

}  // namespace reshala
