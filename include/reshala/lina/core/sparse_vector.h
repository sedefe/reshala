#pragma once

#include <assert.h>

#include <algorithm>
#include <ostream>
#include <vector>

#include "reshala/constants.h"
#include "reshala/lina/core/dense_vector.h"
#include "reshala/numerics.h"
#include "reshala/utils.h"

namespace reshala {

class SparseVector {
   public:
    class IteratorBase;

    SparseVector(Index dim) : dim_(dim) {}
    SparseVector(Index dim, Index i, Scalar v) : dim_(dim) { Push(i, v); }
    SparseVector(Index dim, const std::vector<Index> &indices, const std::vector<Scalar> &values,
                 bool sorted) {
        dim_ = dim;
        indices_ = std::move(indices);
        values_ = std::move(values);
        if (!sorted) {
            Sort();
        }
    }
    SparseVector(Index dim, const Scalar *array) : dim_(dim) {
        for (Index i = 0; i < dim; i++) {
            if (!IsZero(array[i])) {
                Push(i, array[i]);
            }
        }
    }
    SparseVector(const DenseVector &vec) : SparseVector(vec.size(), vec.data()) {}
    DenseVector ToDense() const;

    Scalar At(Index i) const {
        auto offset = FindOffset(i);
        return (offset != Size() && indices_.data()[offset] == i) ? values_[offset] : Scalar(0);
    }
    Scalar &AtRef(Index i) {
        auto offset = FindOffset(i);
        assert(offset != Size());
        return values_[offset];
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

    void EraseOffset(Index offset) {  // Erase by offset
        assert(offset >= 0 && offset < Size());
        indices_.erase(indices_.begin() + offset);
        values_.erase(values_.begin() + offset);
    }

    void EraseIndex(Index i) {  // Erase by index, no reindexing
        auto offset = FindOffset(i);
        EraseOffset(offset);  // Reuse the offset version
    }

    inline Index FindOffset(Index i) const {
        return std::lower_bound(indices_.data(), indices_.data() + indices_.size(), i) -
               indices_.data();
    }

    void Insert(Index i, Scalar v, Index offset) {
        indices_.insert(indices_.begin() + offset, i);
        values_.insert(values_.begin() + offset, v);
    }

    void Sort();

    inline Scalar Norm2() const {
        Scalar res = 0.0;
        for (auto v : values_) res += v * v;
        return res;
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
