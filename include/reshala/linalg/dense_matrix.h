#pragma once

#include <assert.h>

#include <vector>

#include "reshala/linalg/types.h"

namespace reshala {

class DenseMatrix {
   public:
    DenseMatrix(Index m, Index n) : m_(m), n_(n) { data_.resize(m * n); };

    Index GetNRows() const { return m_; }
    Index GetNCols() const { return n_; }

    void resizeAsZero(Index m, Index n) {
        m_ = m;
        n_ = n;
        data_.assign(m_ * n_, 0.0);
    }

    Scalar* rowView(Index i) {
        assert(i < m_);
        return data_.data() + i * n_;
    }

    const Scalar* rowView(Index i) const {
        assert(i < m_);
        return data_.data() + i * n_;
    }

   private:
    Index m_, n_;
    std::vector<Scalar> data_;
};

}  // namespace reshala
