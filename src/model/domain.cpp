#include "reshala/model/domain.h"

#include <cmath>
#include <cstdio>

namespace reshala {

VarType Bounds2Type(const Bounds &bounds) {
    auto l = bounds.le;
    auto u = bounds.ri;
    if ((-kInf < l) && (u == kInf)) {
        return VarType::kLower;
    } else if ((-kInf == l) && (u < kInf)) {
        return VarType::kUpper;
    } else if ((-kInf < l) && (l + kEpsZero < u) && (u < kInf)) {
        return VarType::kBoxed;
    } else if ((-kInf == l) && (u == kInf)) {
        return VarType::kFree;
    } else if ((-kInf < l) && (std::abs(u - l) <= kEpsZero) && (u < kInf)) {
        return VarType::kFixed;
    }
    printf("Strange bounds for variable: %f .. %f\n", l, u);
    return VarType::kUnknown;
}

}  // namespace reshala
