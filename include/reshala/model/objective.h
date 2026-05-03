#pragma once

#include "reshala/linalg/dense_vector.h"
#include "reshala/linalg/operators.h"
#include "reshala/linalg/sparse_vector.h"

namespace reshala {

enum class Sense { kMin, kMax };

class Objective {
   public:
    DenseVector coefficients;  // Coefficients for decision variables
    Scalar c0;                 // Constant term
    Sense sense;               // Minimize or maximize

    Objective() : coefficients(), c0(0.0), sense(Sense::kMin) {}

    Objective(const Objective& other) = default;
    Objective(Objective&& other) noexcept = default;

    Objective& operator=(const Objective& other) = default;
    Objective& operator=(Objective&& other) noexcept = default;

    Scalar evaluate(const DenseVector& x) const {
        Scalar res = 0;
        dot(coefficients, x, res);
        return res + c0;
    }

    Scalar evaluate(const SparseVector& x) const {
        Scalar res = 0;
        dot(coefficients, x, res);
        return res + c0;
    }

    friend std::ostream& operator<<(std::ostream& os, const Objective& obj) {
        os << "Objective: ";
        if (obj.sense == Sense::kMin) {
            os << "Minimize ";
        } else {
            os << "Maximize ";
        }

        // Print coefficients
        os << obj.coefficients << " * x + " << obj.c0;
        return os;
    }
};

}  // namespace reshala
