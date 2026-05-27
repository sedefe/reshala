#pragma once

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

}  // namespace reshala
