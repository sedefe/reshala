#pragma once

#include <array>

#include "reshala/milp/utils.h"
#include "reshala/utils.h"

namespace reshala {

class History {
   public:
    History(Index n_vars) : s_(n_vars, {0.0, 0.0}), s2_(n_vars, {0.0, 0.0}), n_(n_vars, {0, 0}) {}

    void Add(Index iv, Direction dir, Scalar dy, Scalar dx);
    Scalar Estimate(Index iv, Direction dir, Scalar dx) const;
    Scalar GetSigma(Index iv, Direction dir) const;
    inline Index GetN(Index iv, Direction dir) const { return n_[iv][Dir2Index(dir)]; }

   private:
    std::vector<std::array<Scalar, 2>> s_;
    std::vector<std::array<Scalar, 2>> s2_;
    std::vector<std::array<Index, 2>> n_;
};

}  // namespace reshala
