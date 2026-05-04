#pragma once

#include <vector>

#include "reshala/linalg/constants.h"

namespace reshala {

struct Bounds {
    Scalar le = -kInf;
    Scalar ri = kInf;
};

inline Bounds BoundsIntersection(const Bounds &bnd1, const Bounds &bnd2) {
    return Bounds{std::max(bnd1.le, bnd2.le), std::min(bnd1.ri, bnd2.ri)};
}

inline bool InBounds(Scalar val, const Bounds &bounds, Scalar eps) {
    return val + eps >= bounds.le && val - eps <= bounds.ri;
}

enum class VarType {
    kLower,
    kFree,
    kFixed,
    kUpper,
    kBoxed,
    kUnknown
};  // Default is lower, because default bounds are [0; inf)
VarType GetType(const Bounds &bounds);

struct Domain {
    std::vector<Bounds> bounds;
    std::vector<VarType> types;
    std::vector<bool> integrality;

    inline void SetVarBounds(Index iv, const Bounds &bnds) {
        bounds[iv] = bnds;
        types[iv] = GetType(bnds);
    }

    void Resize(Index n) {
        bounds.resize(n);
        types.resize(n);
        integrality.resize(n);
    }

    size_t Size() const { return bounds.size(); }
    void Push(const Bounds &b, bool i) {
        bounds.push_back(b);
        types.push_back(GetType(b));
        integrality.push_back(i);
    }
};

}  // namespace reshala
