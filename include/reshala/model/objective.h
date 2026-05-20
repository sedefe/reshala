#pragma once

#include "reshala/linalg/dense_vector.h"
#include "reshala/linalg/operators.h"
#include "reshala/linalg/sparse_vector.h"

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
        return res * mult + c0;
    }

    Scalar evaluate(const SparseVector& x) const {
        Scalar res = 0;
        dot(coefficients, x, res);
        return res * mult + c0;
    }

    friend std::ostream& operator<<(std::ostream& os, const Objective& obj) {
        Scalar sense = 0;
        bool first_coeff = true;
        if (obj.orig_sense == Sense::kMin) {
            os << "Minimize\n";
            sense = 1.;
        } else {
            os << "Maximize\n";
            sense = -1.;
        }
        for (Index i = 0; i < obj.coefficients.size(); i++) {
            if (obj.coefficients[i] != 0) {
                if (!first_coeff) {
                    os << " + ";
                }
                first_coeff = false;
                os << sense * obj.coefficients[i] << " x" << i;
            }
        }
        if (obj.c0 != 0) {
            os << " + " << obj.c0;
        }

        return os;
    }
};

}  // namespace reshala
