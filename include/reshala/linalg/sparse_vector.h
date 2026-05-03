#pragma once

#include <algorithm>
#include <cassert>
#include <ostream>
#include <vector>

#include "constants.h"

namespace reshala {

class SparseVector {
   public:
    SparseVector(Index dim) : dim_(dim) {}

    // Build from unsorted or sorted lists (caller guarantees sorted & unique? we sort)
    SparseVector(Index dim, const std::vector<Index> &indices, const std::vector<Scalar> &values) {
        dim_ = dim;
        assign(indices, values);
    }

    Scalar get(Index idx) const {
        auto pos = findIndex(idx);
        return (pos != indices_.end()) ? values_[pos - indices_.begin()] : Scalar(0);
    }

    size_t size() const { return indices_.size(); }
    void clear() {
        indices_.clear();
        values_.clear();
    }

    void reserve(size_t n) {
        indices_.reserve(n);
        values_.reserve(n);
    }

    void push(Index i, Scalar x) {
        indices_.push_back(i);
        values_.push_back(x);
    }

    const std::vector<Index> &indices() const { return indices_; }
    const std::vector<Scalar> &values() const { return values_; }
    Index dim() const { return dim_; }

    // ----- modify whole structure -----
    void assign(const std::vector<Index> &indices, const std::vector<Scalar> &values) {
        // copy and sort by index, merging duplicates by summing values
        std::vector<std::pair<Index, Scalar>> pairs;
        pairs.reserve(indices.size());
        for (size_t i = 0; i < indices.size(); ++i) pairs.emplace_back(indices[i], values[i]);
        std::sort(pairs.begin(), pairs.end(),
                  [](const auto &a, const auto &b) { return a.first < b.first; });

        indices_.clear();
        values_.clear();
        for (const auto &p : pairs) {
            if (!indices_.empty() && indices_.back() == p.first) {
                values_.back() += p.second;
                if (values_.back() == Scalar(0)) {
                    indices_.pop_back();
                    values_.pop_back();
                }
            } else if (p.second != Scalar(0)) {
                indices_.push_back(p.first);
                values_.push_back(p.second);
            }
        }
    }

    friend std::ostream &operator<<(std::ostream &os, const SparseVector &sv);

   private:
    Index dim_;
    std::vector<Index> indices_;
    std::vector<Scalar> values_;

    typename std::vector<Index>::const_iterator findIndex(Index idx) const {
        return std::lower_bound(indices_.begin(), indices_.end(), idx);
    }
    typename std::vector<Index>::iterator findIndex(Index idx) {
        return std::lower_bound(indices_.begin(), indices_.end(), idx);
    }

    void erase(size_t pos) {
        indices_.erase(indices_.begin() + pos);
        values_.erase(values_.begin() + pos);
    }

    void insert(Index idx, Scalar val, typename std::vector<Index>::iterator pos) {
        size_t i = pos - indices_.begin();
        indices_.insert(pos, idx);
        values_.insert(values_.begin() + i, val);
    }
};

}  // namespace reshala
