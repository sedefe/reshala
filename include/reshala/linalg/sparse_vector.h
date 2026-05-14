#pragma once

#include <assert.h>

#include <algorithm>
#include <ostream>
#include <vector>

#include "reshala/constants.h"
#include "reshala/utils.h"

namespace reshala {

class SparseVector {
   public:
    class Iterator;

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

    Scalar At(Index idx) const {
        auto pos = FindIndex(idx);
        return (pos != indices_.end()) ? values_[pos - indices_.begin()] : Scalar(0);
    }

    size_t Size() const { return indices_.size(); }
    void Clear() {
        indices_.clear();
        values_.clear();
    }

    void Resize(size_t n) {
        indices_.resize(n);
        values_.resize(n);
    }

    void Reserve(size_t n) {
        indices_.reserve(n);
        values_.reserve(n);
    }

    void Push(Index i, Scalar x) {
        indices_.push_back(i);
        values_.push_back(x);
    }

    void Erase(Index i) {
        auto it = std::lower_bound(indices_.begin(), indices_.end(), i);
        assert(it != indices_.end() && *it == i);

        Index pos = std::distance(indices_.begin(), it);
        indices_.erase(indices_.begin() + pos);
        values_.erase(values_.begin() + pos);

        for (Index j = pos; j < indices_.size(); ++j) {
            indices_[j]--;
        }
    }

    const std::vector<Index> &indices() const { return indices_; }
    std::vector<Index> &indices() { return indices_; }
    const std::vector<Scalar> &values() const { return values_; }
    std::vector<Scalar> &values() { return values_; }
    Index dim() const { return dim_; }
    void SetDim(Index dim) { dim_ = dim; }

    friend std::ostream &operator<<(std::ostream &os, const SparseVector &sv);

   private:
    Index dim_;
    std::vector<Index> indices_;
    std::vector<Scalar> values_;

    typename std::vector<Index>::const_iterator FindIndex(Index idx) const {
        return std::lower_bound(indices_.begin(), indices_.end(), idx);
    }
    typename std::vector<Index>::iterator FindIndex(Index idx) {
        return std::lower_bound(indices_.begin(), indices_.end(), idx);
    }

    void Insert(Index idx, Scalar val, typename std::vector<Index>::iterator pos) {
        size_t i = pos - indices_.begin();
        indices_.insert(pos, idx);
        values_.insert(values_.begin() + i, val);
    }
};

class SparseVector::Iterator {
   public:
    Iterator(const SparseVector &sv) : sv_(sv), pos_(0) {}

    operator bool() const { return pos_ < sv_.indices_.size(); }

    Iterator &operator++() {
        if (pos_ < sv_.indices_.size()) {
            ++pos_;
        }
        return *this;
    }

    std::pair<Index, Scalar> operator*() const { return {sv_.indices_[pos_], sv_.values_[pos_]}; }

    Index index() const { return sv_.indices_[pos_]; }
    Scalar value() const { return sv_.values_[pos_]; }

   private:
    const SparseVector &sv_;
    size_t pos_;
};

using SvIterator = SparseVector::Iterator;

}  // namespace reshala
