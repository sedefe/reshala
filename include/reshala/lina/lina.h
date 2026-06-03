#pragma once

#include "reshala/lina/dense_matrix.h"
#include "reshala/lina/sparse_matrix.h"

namespace reshala {

class Lina {
   public:
    Lina(Index m = 0, Index n = 0) : m_(m), n_(n), Binv_(m, m) {}
    Lina(const Lina& other) : m_(other.m_), n_(other.n_), Binv_(other.Binv_) {}

    void Init();

    void Btran(Index iv, DenseVector& res);
    void Ftran(const SparseVector& sv, DenseVector& res);
    void Update(Index iv, const SparseVector& delta);

   private:
    Index m_;
    Index n_;
    DenseMatrix Binv_;
};

}  // namespace reshala
