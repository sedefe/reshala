#pragma once

#include <array>

#include "reshala/milp/utils.h"
#include "reshala/utils.h"

namespace reshala {

class History {
    Index kMinSamples = 8;

   public:
    History(Index n_vars) : s_(n_vars, {0.0, 0.0}), s2_(n_vars, {0.0, 0.0}), n_(n_vars, {0, 0}) {}

    void Add(Index iv, Direction dir, Scalar dy, Scalar dx);

    inline Scalar Estimate(Index iv, Direction dir, Scalar dx) const {
        return GetMean(iv, dir) * dx;
    }

    inline Scalar GetMean(Index iv, Direction dir) const {
        Index d = Dir2Index(dir);
        return s_[iv][d] / n_[iv][d];
    }

    inline Scalar GetVariance(Index iv, Direction dir) const {
        Index d = Dir2Index(dir);
        return (s2_[iv][d] - s_[iv][d] * s_[iv][d] / n_[iv][d]) / (n_[iv][d] - 1);
    }

    bool IsEnough(Index iv) const;

   private:
    std::vector<std::array<Scalar, 2>> s_;
    std::vector<std::array<Scalar, 2>> s2_;
    std::vector<std::array<Index, 2>> n_;
};

}  // namespace reshala
