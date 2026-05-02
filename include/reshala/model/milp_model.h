#pragma once

#include "reshala/linalg/sparse_matrix.h"
#include "reshala/model/objective.h"

namespace reshala {

struct Bound {
    Scalar le;
    Scalar ri;
};

class MilpModel {
   public:
    MilpModel() = default;

   private:
    Objective obj;
    SparseColMatrix Ac;
    SparseRowMatrix Ar;
    std::vector<Bound> rhs;
    std::vector<Bound> bounds;
};

}  // namespace reshala
