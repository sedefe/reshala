#include "reshala/lina/sparse_vector.h"

namespace reshala {

template <typename Operation>
SparseVector combine(const SparseVector& sv1, const SparseVector& sv2, Operation op) {
    assert(sv1.dim() == sv2.dim() && "SparseVector combine: vectors are of different dimensions");

    std::vector<Index> result_indices;
    std::vector<Scalar> result_values;

    const auto& x_indices = sv1.indices();
    const auto& x_values = sv1.values();
    const auto& y_indices = sv2.indices();
    const auto& y_values = sv2.values();

    size_t i = 0, j = 0;
    size_t x_size = x_indices.size();
    size_t y_size = y_indices.size();

    while (i < x_size || j < y_size) {
        Index idx;
        Scalar x_val = Scalar(0);
        Scalar y_val = Scalar(0);

        if (i < x_size && (j >= y_size || x_indices[i] < y_indices[j])) {
            idx = x_indices[i];
            x_val = x_values[i];
            i++;
        } else if (j < y_size && (i >= x_size || y_indices[j] < x_indices[i])) {
            idx = y_indices[j];
            y_val = y_values[j];
            j++;
        } else {
            idx = x_indices[i];
            x_val = x_values[i];
            y_val = y_values[j];
            i++;
            j++;
        }

        Scalar result_val = op(x_val, y_val);

        if (!IsZero(result_val)) {
            result_indices.push_back(idx);
            result_values.push_back(result_val);
        }
    }

    return SparseVector(sv1.dim(), result_indices, result_values);
}

SparseVector operator+(const SparseVector& sv1, const SparseVector& sv2) {
    return combine(sv1, sv2, [](Scalar x, Scalar y) { return x + y; });
}

SparseVector operator-(const SparseVector& sv1, const SparseVector& sv2) {
    return combine(sv1, sv2, [](Scalar x, Scalar y) { return x - y; });
}

SparseVector operator*(SparseVector sv, Scalar x) {
    sv *= x;
    return sv;
}

SparseVector operator*(Scalar x, SparseVector sv) {
    sv *= x;
    return sv;
}

std::ostream& operator<<(std::ostream& os, const SparseVector& sv) {
    if (sv.indices_.empty()) {
        os << "0";
        return os;
    }

    for (Index i = 0; i < sv.indices_.size(); ++i) {
        Scalar v = sv.values_[i];
        if (i == 0) {
            os << (v >= 0 ? "" : "- ");
        } else {
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
