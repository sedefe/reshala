#include "reshala/linalg/dense_vector.h"

namespace reshala {

std::ostream &operator<<(std::ostream &os, const DenseVector &dv) {
    bool first = true;
    for (size_t i = 0; i < dv.size(); ++i) {
        if (!IsZero(dv[i])) {
            if (!first) {
                os << " + ";
            }
            os << dv[i] << " * x[" << i << "]";
            first = false;
        }
    }
    if (first) {
        os << "0";
    }
    os << "\n";

    return os;
}

}  // namespace reshala
