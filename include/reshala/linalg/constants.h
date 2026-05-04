#pragma once

#include <limits>
#include "types.h"

namespace reshala {

constexpr Scalar kEpsZero = 1e-8;
constexpr Scalar kInf = std::numeric_limits<Scalar>::max();
constexpr Scalar kMaxAbs = 1e6;


inline bool isZero(Scalar x) {
    return x < kEpsZero and x > -kEpsZero;
}

};
