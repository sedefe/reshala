#pragma once

#include <limits>

#include "reshala/types.h"

namespace reshala {

constexpr Index kMaxInt = std::numeric_limits<Index>::max();

constexpr Scalar kEpsZero = 1e-6;
constexpr Scalar kInf = std::numeric_limits<Scalar>::infinity();
constexpr Scalar kNan = std::numeric_limits<Scalar>::quiet_NaN();
constexpr Scalar kMaxAbs = 1e6;
constexpr Scalar kMipGap = 1e-4;

};  // namespace reshala
