#pragma once

#include <assert.h>

#include <algorithm>
#include <ostream>
#include <vector>

#include "reshala/constants.h"
#include "reshala/linalg/dense_vector.h"
#include "reshala/utils.h"

namespace reshala {

class SparseVector {
   public:
    class IteratorBase;

    SparseVector(Index dim) : dim_(dim) {}
    SparseVector(Index dim, Index i, Scalar v) : dim_(dim) { Push(i, v); }
    SparseVector(Index dim, const std::vector<Index> &indices, const std::vector<Scalar> &values) {
        dim_ = dim;
        indices_ = std::move(indices);
        values_ = std::move(values);
    }
    SparseVector(Index dim, const Scalar *array) : dim_(dim) {
        for (Index i = 0; i < dim; i++) {
            if (!IsZero(array[i])) {
                Push(i, array[i]);
            }
        }
    }
    SparseVector(const DenseVector &vec) : SparseVector(vec.size(), vec.data()) {}

    // Scalar At(Index i) const {
    //     auto pos = FindIndex(i);
    //     return (pos != indices_.end()) ? values_[pos - indices_.begin()] : Scalar(0);
    // }
    Scalar &AtRef(Index i) {
        auto pos = FindIndex(i);
        assert(pos != indices_.end());
        return values_[pos - indices_.begin()];
    }

    inline size_t Size() const { return indices_.size(); }
    inline void Clear() {
        indices_.clear();
        values_.clear();
    }

    inline bool Empty() const { return indices_.empty(); }

    inline void Resize(size_t n) {
        indices_.resize(n);
        values_.resize(n);
    }

    inline void Reserve(size_t n) {
        indices_.reserve(n);
        values_.reserve(n);
    }

    inline void Push(Index i, Scalar x) {
        indices_.push_back(i);
        values_.push_back(x);
    }

    void Erase(typename std::vector<Index>::iterator it) {  // Erase by iterator
        // Verify the iterator is valid (optional, for debugging)
        assert(it >= indices_.begin() && it < indices_.end());

        Index pos = std::distance(indices_.begin(), it);
        indices_.erase(it);
        values_.erase(values_.begin() + pos);
    }

    void Erase(Index i) {  // Erase by index, no reindexing
        auto it = std::lower_bound(indices_.begin(), indices_.end(), i);
        assert(it != indices_.end() && *it == i);
        Erase(it);  // Reuse the iterator version
    }

    typename std::vector<Index>::iterator FindIndex(Index i) {
        return std::lower_bound(indices_.begin(), indices_.end(), i);
    }

    void Insert(Index i, Scalar v, typename std::vector<Index>::const_iterator pos) {
        Index d = pos - indices_.begin();
        indices_.insert(pos, i);
        values_.insert(values_.begin() + d, v);
    }

    const std::vector<Index> &indices() const { return indices_; }
    std::vector<Index> &indices() { return indices_; }
    const std::vector<Scalar> &values() const { return values_; }
    std::vector<Scalar> &values() { return values_; }
    Index dim() const { return dim_; }
    void SetDim(Index dim) { dim_ = dim; }

    // Arithmetics
    inline SparseVector &operator*=(Scalar x) {
        for (Scalar &v : values_) v *= x;
        return *this;
    }
    friend std::ostream &operator<<(std::ostream &os, const SparseVector &sv);

   private:
    Index dim_;
    std::vector<Index> indices_;
    std::vector<Scalar> values_;
};

SparseVector operator+(const SparseVector &sv1, const SparseVector &sv2);
SparseVector operator-(const SparseVector &sv1, const SparseVector &sv2);
SparseVector operator*(SparseVector sv, Scalar x);
SparseVector operator*(Scalar x, SparseVector sv);

template <bool IsConst>
class IteratorBase {
    using VectorType = std::conditional_t<IsConst, const SparseVector, SparseVector>;
    using ScalarType = std::conditional_t<IsConst, const Scalar, Scalar>;
    using IndexType = std::conditional_t<IsConst, const Index, Index>;

   public:
    IteratorBase(VectorType &sv) : sv_(sv), pos_(0) {}

    operator bool() const { return pos_ < sv_.indices().size(); }

    IteratorBase &operator++() {
        if (pos_ < sv_.indices().size()) {
            ++pos_;
        }
        return *this;
    }

    std::pair<Index, Scalar> operator*() const { return {sv_.indices()[pos_], sv_.values()[pos_]}; }

    inline Index index() const { return sv_.indices()[pos_]; }
    inline IndexType &indexRef() { return sv_.indices()[pos_]; }
    inline Scalar value() const { return sv_.values()[pos_]; }
    inline ScalarType &valueRef() { return sv_.values()[pos_]; }

   private:
    VectorType &sv_;
    Index pos_;
};

using MutableSvIterator = IteratorBase<false>;
using SvIterator = IteratorBase<true>;

}  // namespace reshala
