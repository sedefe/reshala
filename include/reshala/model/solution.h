#pragma once

#include <string>
#include <vector>

#include "reshala/constants.h"

namespace reshala {

enum class LpStatus { kUnknown, kOptimal, kInfeasible, kUnbounded };
std::string LpStatus2Str(LpStatus status);

struct Solution {
    LpStatus status = LpStatus::kUnknown;
    Scalar y = kInf;
    std::vector<Scalar> x;
};

inline Solution InfeasibleSolution() { return {LpStatus::kInfeasible, kInf, {}}; }

}  // namespace reshala
