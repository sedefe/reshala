#include "reshala/linalg/sparse_vector.h"

namespace reshala {

std::ostream& operator<<(std::ostream& os, const SparseVector& sv) {
    if (sv.indices_.empty()) {
        os << "0";
        return os;
    }

    for (size_t i = 0; i < sv.indices_.size(); ++i) {
        if (i > 0) {
            os << " + ";
        }
        os << sv.values_[i] << " x" << sv.indices_[i];
    }

    return os;
}

};  // namespace reshala
