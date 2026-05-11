#pragma once

#include <chrono>
#include <cmath>

#include "reshala/constants.h"

namespace reshala {

inline Scalar GetFraction(Scalar x) {
    Scalar nearest = std::round(x);
    return std::abs(x - nearest);
}

}  // namespace reshala
