#include "reshala/lina/core/sparse_vector.h"

#include <numeric>

namespace reshala {

DenseVector SparseVector::ToDense() const {
    DenseVector res(dim_, 0.0);
    for (SvIterator el(*this); el; ++el) {
        res[el.index()] = el.value();
    }
    return res;
}

void SparseVector::Sort() {
    Index n = indices_.size();
    if (n <= 1) return;

    // Create permutation indices
    std::vector<size_t> perm(n);
    std::iota(perm.begin(), perm.end(), 0);

    // Sort permutation based on indices
    std::sort(perm.begin(), perm.end(),
              [this](Index i, Index j) { return indices_[i] < indices_[j]; });

    // Apply permutation to both vectors
    std::vector<Index> sorted_indices(n);
    std::vector<Scalar> sorted_values(n);

    for (Index i = 0; i < n; ++i) {
        sorted_indices[i] = indices_[perm[i]];
        sorted_values[i] = values_[perm[i]];
    }

    indices_.swap(sorted_indices);
    values_.swap(sorted_values);
}

template <typename Operation>
SparseVector combine(const SparseVector& sv1, const SparseVector& sv2, Operation op) {
    assert(sv1.dim() == sv2.dim() && "SparseVector combine: vectors are of different dimensions");

    SparseVector res(sv1.dim());

    const auto& ind1 = sv1.indices();
    const auto& val1 = sv1.values();
    const auto& ind2 = sv2.indices();
    const auto& val2 = sv2.values();

    Index i1 = 0, i2 = 0;
    auto n1 = sv1.Size();
    auto n2 = sv2.Size();
    res.Reserve(n1 + n2);

    while (i1 < n1 && i2 < n2) {
        Index ind;
        Scalar v1 = Scalar(0);
        Scalar v2 = Scalar(0);

        if (ind1[i1] < ind2[i2]) {
            ind = ind1[i1];
            v1 = val1[i1];
            i1++;
        } else if (ind2[i2] < ind1[i1]) {
            ind = ind2[i2];
            v2 = val2[i2];
            i2++;
        } else {
            ind = ind1[i1];
            v1 = val1[i1];
            v2 = val2[i2];
            i1++;
            i2++;
        }

        Scalar val = op(v1, v2);
        if (!IsZero(val)) {
            res.Push(ind, val);
        }
    }
    for (; i1 < n1; i1++) {
        res.Push(ind1[i1], val1[i1]);
    }
    for (; i2 < n2; i2++) {
        res.Push(ind2[i2], op(0, val2[i2]));
    }

    return res;
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
