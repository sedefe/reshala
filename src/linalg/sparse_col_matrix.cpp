#include <assert.h>

#include "reshala/linalg/sparse_matrix.h"

namespace reshala {

void MulScmSv(const SparseColMatrix &scm, const SparseVector &sv, DenseVector &res) {
    auto m = scm.GetNRows();
    auto n = scm.GetNCols();
    assert(n == sv.dim());
    res.assign(m, Scalar(0));

    for (Index i = 0; i < sv.size(); i++) {
        auto ind = sv.indices()[i];
        auto val = sv.values()[i];
        const auto &col = scm.getCols()[ind];
        for (Index j = 0; j < col.size(); j++) {
            res[col.indices()[j]] += val * col.values()[j];
        }
    }
}

}  // namespace reshala
