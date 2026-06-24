#pragma once

#include <chrono>

#include "reshala/constants.h"
#include "reshala/logging.h"

namespace reshala {

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

#if defined(NDEBUG)
    const std::string kBuildType = "release";
#elif defined(DEBUG)
    const std::string kBuildType = "debug";
#else
    const std::string kBuildType = "X3";
#endif

}  // namespace reshala
