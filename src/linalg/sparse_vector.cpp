#include "reshala/linalg/sparse_vector.h"

namespace reshala {

std::ostream& operator<<(std::ostream& os, const SparseVector& sv) {
    if (sv.indices_.empty()) {
        os << "0";
        return os;
    }

    for (Index i = 0; i < sv.indices_.size(); ++i) {
        Scalar v = sv.values_[i];
        if (i > 0) {
            os << (v >= 0 ? " + " : " - ");
        }
        v = std::abs(v);
        if (v == 1) {
            os << "x" << sv.indices_[i];
        } else {
            os << v << " x" << sv.indices_[i];
        }
    }

    return os;
}

};  // namespace reshala
