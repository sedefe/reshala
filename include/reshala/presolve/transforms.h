#pragma once

#include "reshala/lina/core/sparse_vector.h"
#include "reshala/model/solution.h"

namespace reshala {

class Transform {
   public:
    virtual ~Transform() = default;

    virtual void Undo(Solution& sol) = 0;
};

class FixVariableTransform : public Transform {
   public:
    FixVariableTransform(Index iv, Scalar val) : iv_(iv), val_(val) {}

    void Undo(Solution& sol) override;

   private:
    Index iv_;
    Scalar val_;
};

class SimpleSubTransform : public Transform {  // iv1 <- a * iv2 + b
   public:
    SimpleSubTransform(Index iv1, Scalar a, Index iv2, Scalar b)
        : iv1_(iv1), a_(a), iv2_(iv2), b_(b) {}

    void Undo(Solution& sol) override;

   private:
    Index iv1_;
    Scalar a_;
    Index iv2_;
    Scalar b_;
};

class LinCombTransform : public Transform {  // iv <- Sum(a_k * x_k) + b
   public:
    LinCombTransform(Index iv, const SparseVector& sv, Scalar b) : iv_(iv), sv_(sv), b_(b) {}

    void Undo(Solution& sol) override;

   private:
    Index iv_;
    SparseVector sv_;
    Scalar b_;
};

}  // namespace reshala
