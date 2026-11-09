#pragma once

#include "reshala/lina/core/operators.h"

namespace reshala {

enum class Sense { kMin, kMax };

class Objective {
   public:
    DenseVector coefficients;
    Scalar mult;
    Scalar c0;
    Sense orig_sense;  // for output only, we actually always minimize

    Objective() : coefficients(), mult(1.0), c0(0.0), orig_sense(Sense::kMin) {}

    Objective(const Objective& other) = default;
    Objective(Objective&& other) noexcept = default;

    Objective& operator=(const Objective& other) = default;
    Objective& operator=(Objective&& other) noexcept = default;

    Scalar evaluate(const DenseVector& x) const {
        Scalar res = 0;
        dot(coefficients, x, res);
        return (res + c0) * mult;
    }

    Scalar evaluate(const SparseVector& x) const {
        Scalar res = 0;
        dot(coefficients, x, res);
        return (res + c0) * mult;
    }
};

}  // namespace reshala
