#pragma once

#include <vector>

#include "reshala/constants.h"

namespace reshala {

struct Bounds {
    Scalar le;
    Scalar ri;

    Bounds(Scalar l = 0.0, Scalar r = kInf) : le(l), ri(r) {}
};

inline Bounds BoundsIntersection(const Bounds &bnd1, const Bounds &bnd2) {
    return Bounds{std::max(bnd1.le, bnd2.le), std::min(bnd1.ri, bnd2.ri)};
}

inline bool InBounds(Scalar val, const Bounds &bounds, Scalar eps) {
    return val + eps >= bounds.le && val - eps <= bounds.ri;
}

enum class BndType { kLower, kFree, kFixed, kUpper, kBoxed, kInfeasible };
BndType Bounds2Type(const Bounds &bounds);

class Domain {
   public:
    inline const std::vector<Bounds> &GetBounds() const { return bounds_; }
    inline const Bounds &GetBounds(Index iv) const { return bounds_[iv]; }
    inline void SetBounds(Index iv, const Bounds &bnds) {
        bounds_[iv] = bnds;
        types_[iv] = Bounds2Type(bnds);
    }
    inline const BndType &GetType(Index iv) const { return types_[iv]; }

    inline bool GetIntegrality(Index iv) const { return integrality_[iv]; }
    inline void SetIntegrality(Index iv, bool b) { integrality_[iv] = b; }

    inline size_t Size() const { return bounds_.size(); }
    inline void Resize(Index n) {
        bounds_.resize(n);
        types_.resize(n);
        integrality_.resize(n);
    }
    inline void Push(const Bounds &b, bool i) {
        bounds_.push_back(b);
        types_.push_back(Bounds2Type(b));
        integrality_.push_back(i);
    }
    inline void Move(Index i_read, Index i_write) {
        bounds_[i_write] = std::move(bounds_[i_read]);
        types_[i_write] = std::move(types_[i_read]);
        integrality_[i_write] = std::move(integrality_[i_read]);
    }

   private:
    std::vector<Bounds> bounds_;
    std::vector<BndType> types_;
    std::vector<bool> integrality_;
};

}  // namespace reshala
