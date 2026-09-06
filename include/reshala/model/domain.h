#pragma once

#include <array>
#include <vector>

#include "reshala/constants.h"
#include "reshala/lina/core/masked_vector.h"

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

    inline bool GetIntegrality(Index iv) const { return integrality_.Get(iv); }
    inline void SetIntegrality(Index iv, bool b) { integrality_.Set(iv, b); }

    inline size_t Size() const { return bounds_.size(); }
    inline void Resize(Index n) {
        bounds_.resize(n);
        types_.resize(n);
        integrality_.Resize(n);
    }
    inline void Push(const Bounds &bnd, bool is_int) {
        bounds_.push_back(bnd);
        types_.push_back(Bounds2Type(bnd));
        integrality_.Push(is_int);
    }
    inline void Move(Index i_read, Index i_write) {
        bounds_[i_write] = std::move(bounds_[i_read]);
        types_[i_write] = std::move(types_[i_read]);
        integrality_.Set(i_write, integrality_.Get(i_read));
    }

   private:
    std::vector<Bounds> bounds_;
    std::vector<BndType> types_;
    BitMask integrality_;
};

enum class LockType { kDown = 0, kUp = 1 };
inline Index LockType2Index(LockType lt) { return static_cast<Index>(lt); }

struct Locks {
    std::vector<std::array<Scalar, 2>> n_locks;

    inline void Resize(Index n) { n_locks.resize(n, {0, 0}); }
};

}  // namespace reshala
