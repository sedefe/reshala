#pragma once

#include <cmath>
#include <numeric>
#include <vector>

#include "reshala/constants.h"

namespace reshala {

inline bool IsZero(Scalar x) { return x < kEpsZero and x > -kEpsZero; }

inline bool WeakGe(Scalar x, Scalar y) { return x - y > -kEpsZero; }
inline bool StrongGt(Scalar x, Scalar y) { return x - y > kEpsZero; }
inline bool WeakLe(Scalar x, Scalar y) { return x - y < kEpsZero; }
inline bool StrongLt(Scalar x, Scalar y) { return x - y < -kEpsZero; }
inline bool WeakEq(Scalar x, Scalar y) { return IsZero(x - y); }

inline Scalar Floor(Scalar x) { return std::floor(x); }
inline Scalar WeakFloor(Scalar x) { return std::floor(x + kEpsZero); }
inline Scalar Round(Scalar x) { return std::round(x); }
inline Scalar Ceil(Scalar x) { return std::ceil(x); }
inline Scalar WeakCeil(Scalar x) { return std::ceil(x - kEpsZero); }

inline Scalar Fraction(Scalar x) { return x - Floor(x); }

inline Scalar MinFraction(Scalar x) {
    Scalar nearest = std::round(x);
    return std::abs(x - nearest);
}

inline Index GetGcd(const std::vector<Scalar>& vec) {
    Index gcd = 0;
    for (auto x : vec) {
        if ((std::abs(x) <= kMaxInt) and MinFraction(x) == 0.0) {
            gcd = std::gcd(gcd, Index(x));
        } else {
            return 0;
        }
    }
    return gcd;
}

inline void GetExpRange(const std::vector<Scalar>& vec, Index& min_exp, Index& max_exp) {
    min_exp = std::numeric_limits<Index>::max();
    max_exp = std::numeric_limits<Index>::min();
    Index exp;

    for (auto x : vec) {
        if (x == 0) continue;
        std::frexp(x, &exp);
        min_exp = std::min(min_exp, exp);
        max_exp = std::max(max_exp, exp);
    }
}

}  // namespace reshala
