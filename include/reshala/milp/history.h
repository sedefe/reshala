#pragma once

#include <array>

#include "reshala/milp/utils.h"
#include "reshala/utils.h"

namespace reshala {

class History {
   public:
    History(Index n_vars) : s_(n_vars), s2_(n_vars), n_(n_vars) {}

    void Add(Direction dir, Index iv, Scalar dy, Scalar dx);

   private:
    std::vector<std::array<Scalar, 2>> s_;
    std::vector<std::array<Scalar, 2>> s2_;
    std::vector<std::array<Index, 2>> n_;
};

}  // namespace reshala
