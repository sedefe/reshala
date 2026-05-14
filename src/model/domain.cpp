#include "reshala/model/domain.h"
#include "reshala/utils.h"

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
    } else if ((-kInf < l) && IsZero(u - l) && (u < kInf)) {
        return VarType::kFixed;
    }
    std::cerr << "Strange bounds for variable: " << l << " .. " << u << "\n";
    return VarType::kUnknown;
}

}  // namespace reshala
