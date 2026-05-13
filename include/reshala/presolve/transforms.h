#pragma once

#include "reshala/presolve/model_info.h"

namespace reshala {

class Transform {
   public:
    virtual ~Transform() = default;

    virtual void Do(ModelInfo& info) = 0;
    virtual void Undo(Solution& sol) = 0;
};

class FixVariableTransform : public Transform {
   public:
    FixVariableTransform(Index iv, Scalar val)
        : iv_(iv), val_(val) {}

    void Do(ModelInfo& info) override;
    void Undo(Solution& sol) override;

   private:
    Index iv_;
    Scalar val_;
};

}  // namespace reshala
