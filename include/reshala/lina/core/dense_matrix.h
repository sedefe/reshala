#pragma once

#include <assert.h>

#include <vector>

#include "reshala/types.h"

namespace reshala {

class DenseMatrix {
   public:
    DenseMatrix() {}
    DenseMatrix(Index m, Index n) : m_(m), n_(n) { data_.resize(m * n); };

    inline Index GetNRows() const { return m_; }
    inline Index GetNCols() const { return n_; }

    inline void ResizeAsZero(Index m, Index n) {
        m_ = m;
        n_ = n;
        data_.assign(m_ * n_, 0.0);
    }

    inline Scalar* operator[](Index i) {
        assert(i < m_);
        return data_.data() + i * n_;
    }

    inline const Scalar* operator[](Index i) const {
        assert(i < m_);
        return data_.data() + i * n_;
    }

   private:
    Index m_, n_;
    std::vector<Scalar> data_;
};

}  // namespace reshala
