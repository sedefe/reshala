#pragma once

#include <limits>

#include "reshala/types.h"

namespace reshala {

constexpr Scalar kEpsZero = 1e-6;
constexpr Scalar kInf = std::numeric_limits<Scalar>::infinity();
constexpr Scalar kNan = std::numeric_limits<Scalar>::quiet_NaN();
constexpr Scalar kMaxAbs = 1e6;
constexpr Scalar kMipGap = 1e-4;

inline bool IsZero(Scalar x) { return x < kEpsZero and x > -kEpsZero; }

inline bool WeakGe(Scalar x, Scalar y) { return x - y > -kEpsZero; }
inline bool StrongGt(Scalar x, Scalar y) { return x - y > kEpsZero; }
inline bool WeakLe(Scalar x, Scalar y) { return x - y < kEpsZero; }
inline bool StrongLt(Scalar x, Scalar y) { return x - y < -kEpsZero; }
inline bool WeakEq(Scalar x, Scalar y) { return IsZero(x - y); }

};  // namespace reshala
