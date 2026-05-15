#pragma once

namespace reshala {

enum class RuleType { kTrivial, kFast, kMedium, kExhaustive, kUnknown };
enum class RuleResult { kSkipped, kUnchanged, kReduced, kInfeasible, kUnknown };

}  // namespace reshala
