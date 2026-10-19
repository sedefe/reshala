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

inline bool InBounds(Scalar val, const Bounds &bounds) {
    return val + kEpsZero >= bounds.le && val - kEpsZero <= bounds.ri;
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

class Locks {
   public:
    inline void Resize(Index n) {
        n_locks_.resize(n, {0, 0});
        sum_locks_.resize(n, 0);
    }

    inline void Set(Index iv, Index n_down, Index n_up) {
        n_locks_[iv][0] = n_down;
        n_locks_[iv][1] = n_up;
        sum_locks_[iv] = n_down + n_up;
    }
    inline void Inc(Index iv, LockType lt) {
        n_locks_[iv][LockType2Index(lt)]++;
        sum_locks_[iv]++;
    }
    inline void Dec(Index iv, LockType lt) {
        n_locks_[iv][LockType2Index(lt)]--;
        sum_locks_[iv]--;
    }

    inline Index GetNLocks(Index iv, LockType lt) const { return n_locks_[iv][LockType2Index(lt)]; }
    inline Index GetSumLocks(Index iv) const { return sum_locks_[iv]; }

    inline void Move(Index i_read, Index i_write) {
        n_locks_[i_write] = std::move(n_locks_[i_read]);
        sum_locks_[i_write] = std::move(sum_locks_[i_read]);
    }

   private:
    std::vector<std::array<Index, 2>> n_locks_;
    std::vector<Index> sum_locks_;
};

}  // namespace reshala
