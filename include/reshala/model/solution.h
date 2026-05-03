#pragma once

#include <vector>

#include "reshala/linalg/types.h"

namespace reshala {

enum class LpStatus { kOptimal, kInfeasible, kUnbounded };

struct Solution {
    LpStatus status;
    Scalar y;
    std::vector<Scalar> x;
};

}  // namespace reshala
