#pragma once

#include <chrono>
#include <cmath>

#include "reshala/constants.h"
#include "reshala/logging.h"

namespace reshala {

inline bool IsZero(Scalar x) { return x < kEpsZero and x > -kEpsZero; }

inline bool WeakGe(Scalar x, Scalar y) { return x - y > -kEpsZero; }
inline bool StrongGt(Scalar x, Scalar y) { return x - y > kEpsZero; }
inline bool WeakLe(Scalar x, Scalar y) { return x - y < kEpsZero; }
inline bool StrongLt(Scalar x, Scalar y) { return x - y < -kEpsZero; }
inline bool WeakEq(Scalar x, Scalar y) { return IsZero(x - y); }

// Todo test soft floring/ceiling in 3.2/7.2 etc.
inline Scalar Floor(Scalar x) { return std::floor(x); }
inline Scalar WeakFloor(Scalar x) { return std::floor(x + kEpsZero); }
inline Scalar Round(Scalar x) { return std::round(x); }
inline Scalar Ceil(Scalar x) { return std::ceil(x); }
inline Scalar WeakCeil(Scalar x) { return std::ceil(x - kEpsZero); }

inline Scalar GetFraction(Scalar x) {
    Scalar nearest = std::round(x);
    return std::abs(x - nearest);
}

#define MEASURE_TIME(expr)                                                                  \
    [&]() {                                                                                 \
        auto start = std::chrono::high_resolution_clock::now();                             \
        auto result = (expr);                                                               \
        auto end = std::chrono::high_resolution_clock::now();                               \
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start); \
        return std::pair{result, duration.count() / 1e3};                                   \
    }()

static std::string FormatInteger(Index num) {
    if (num == 0) return "0";

    std::string prefix = num < 0 ? "-" : "";
    num = std::abs(num);

    Index scale;
    std::string suffix;
    if (num < 1'000'000) {
        return prefix + std::to_string(num);
    } else if (num < 1'000'000'000) {
        scale = 1000;
        suffix = "k";
    } else {
        scale = 1'000'000;
        suffix = "m";
    }

    return prefix + std::to_string(num / scale) + suffix;
}

}  // namespace reshala
