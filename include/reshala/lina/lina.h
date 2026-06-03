#pragma once

#include "reshala/lina/basis.h"
#include "reshala/lina/core/dense_matrix.h"
#include "reshala/lina/core/sparse_matrix.h"

namespace reshala {

class Lina {
   public:
    Lina() {}
    Lina(const SparseColMatrix* Ac, Index m, Index n, const LpBasis* basis)
        : Ac_(Ac), m_(m), n_(n), basis_(basis), Binv_(m, m) {
        InitBinv();
    }
    Lina& operator=(const Lina&) = default;

    void Btran(Index iv, DenseVector& res);
    void Ftran(Index iv, DenseVector& res);
    void Update(Index iv_leaving, Index iv_entering);

   private:
    const SparseColMatrix* Ac_;
    Index m_;
    Index n_;

    const LpBasis* basis_;
    DenseMatrix Binv_;

    void InitBinv();
};

}  // namespace reshala
