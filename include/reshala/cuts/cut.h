#pragma once

#include "reshala/lina/core/operators.h"

namespace reshala {

struct Cut {
    // lhs >= rhs

    Cut(const SparseVector& l, Scalar r) : lhs(l), rhs(r) {}

    Scalar GetViol(const DenseVector& x) {
        Scalar l;
        dot(lhs, x, l);
        return rhs - l;
    }

    SparseVector lhs;
    Scalar rhs;

    bool selected = false;
    Index age = 0;
};

std::ostream& operator<<(std::ostream& os, const Cut& cut);

}  // namespace reshala
