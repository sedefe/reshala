#pragma once

#include "types.h"

namespace reshala {

constexpr Scalar kEpsZero = 1e-8;

inline bool isZero(Scalar x) {
    return x < kEpsZero and x > -kEpsZero;
}

};
