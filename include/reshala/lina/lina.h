#pragma once

#include "reshala/lina/dense_matrix.h"
#include "reshala/lina/sparse_matrix.h"

namespace reshala {

class Lina {
   public:
    Lina(const SparseColMatrix* Ac = nullptr, Index m = 0, Index n = 0)
        : Ac_(Ac), m_(m), n_(n), Binv_(m, m) {}
    Lina& operator=(const Lina&) = default;

    void Init();
    void SetBasis(const std::vector<Index>&);

    void Btran(Index iv, DenseVector& res);
    void Ftran(Index iv, DenseVector& res);
    void Update(Index iv_leaving, Index iv_entering);

   private:
    const SparseColMatrix* Ac_;
    Index m_;
    Index n_;
    std::vector<Index> basis;
    DenseMatrix Binv_;
};

}  // namespace reshala
