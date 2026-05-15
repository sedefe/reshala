#include "reshala/model/domain.h"
#include "reshala/utils.h"

namespace reshala {

BndType Bounds2Type(const Bounds &bounds) {
    auto l = bounds.le;
    auto u = bounds.ri;
    if ((-kInf < l) && (u == kInf)) {
        return BndType::kLower;
    } else if ((-kInf == l) && (u < kInf)) {
        return BndType::kUpper;
    } else if ((-kInf < l) && (l + kEpsZero < u) && (u < kInf)) {
        return BndType::kBoxed;
    } else if ((-kInf == l) && (u == kInf)) {
        return BndType::kFree;
    } else if ((-kInf < l) && IsZero(u - l) && (u < kInf)) {
        return BndType::kFixed;
    }
    // std::cerr << "Strange bounds for variable: " << l << " .. " << u << "\n";
    return BndType::kInfeasible;
}

}  // namespace reshala
