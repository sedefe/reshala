#pragma once

#include <cmath>
#include <chrono>

#include "reshala/constants.h"
#include "reshala/logging.h"

namespace reshala {

inline bool IsZero(Scalar x) { return x < kEpsZero and x > -kEpsZero; }

inline bool WeakGe(Scalar x, Scalar y) { return x - y > -kEpsZero; }
inline bool StrongGt(Scalar x, Scalar y) { return x - y > kEpsZero; }
inline bool WeakLe(Scalar x, Scalar y) { return x - y < kEpsZero; }
inline bool StrongLt(Scalar x, Scalar y) { return x - y < -kEpsZero; }
inline bool WeakEq(Scalar x, Scalar y) { return IsZero(x - y); }

inline Scalar Floor(Scalar x) { return std::floor(x); }
inline Scalar Round(Scalar x) { return std::round(x); }
inline Scalar Ceil(Scalar x) { return std::ceil(x); }

inline Scalar GetFraction(Scalar x) {
    Scalar nearest = std::round(x);
    return std::abs(x - nearest);
}

}  // namespace reshala
