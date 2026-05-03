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
    MilpModel() {}

    const Objective& getObj() const { return obj; }
    Objective& getObj() { return obj; }

    const SparseColMatrix& getAc() const { return Ac; }
    SparseColMatrix& getAc() { return Ac; }

    const SparseRowMatrix& getAr() const { return Ar; }
    SparseRowMatrix& getAr() { return Ar; }

    const std::vector<Bound>& getRhs() const { return rhs; }
    std::vector<Bound>& getRhs() { return rhs; }

    const std::vector<Bound>& getBounds() const { return bounds; }
    std::vector<Bound>& getBounds() { return bounds; }

   private:
    Objective obj;
    SparseColMatrix Ac;
    SparseRowMatrix Ar;
    std::vector<Bound> rhs;
    std::vector<Bound> bounds;
};

}  // namespace reshala
