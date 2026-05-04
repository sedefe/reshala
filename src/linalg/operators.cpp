#include "reshala/linalg/operators.h"

#include <assert.h>

namespace reshala {

void dot(const SparseVector& sv1, const SparseVector& sv2, Scalar& res) {
    assert(sv1.Size() == sv2.dim());
    res = 0;
}

void dot(const DenseVector& dv1, const SparseVector& sv2, Scalar& res) {
    assert(dv1.size() == sv2.dim());
    res = 0;
    const auto& indices = sv2.indices();
    const auto& values = sv2.values();

    for (size_t i = 0; i < indices.size(); ++i) {
        res += dv1[indices[i]] * values[i];
    }
}

void dot(const DenseVector& dv1, const DenseVector& dv2, Scalar& res) {
    assert(dv1.size() == dv2.size());
    res = 0;
    for (size_t i = 0; i < dv1.size(); ++i) {
        res += dv1[i] * dv2[i];
    }
}

void MulScmSv(const SparseColMatrix& scm, const SparseVector& sv, DenseVector& res) {
    auto m = scm.getNRows();
    auto n = scm.getNCols();
    assert(n == sv.dim());
    res.assign(m, Scalar(0));

    for (Index i = 0; i < sv.Size(); i++) {
        auto ind = sv.indices()[i];
        auto val = sv.values()[i];
        const auto& col = scm.GetCol(ind);
        for (Index j = 0; j < col.Size(); j++) {
            res[col.indices()[j]] += val * col.values()[j];
        }
    }
}

void MulDvDm(const DenseVector& dv, const DenseMatrix& dm, DenseVector& res)  {
    auto m = dm.GetNRows();
    auto n = dm.GetNCols();

    res.assign(n, 0.0);
    for (Index i = 0; i<m; i++) {
        const Scalar * row = dm.RowView(i);
        for (Index j=0; j<n; j++) {
            res[j] += row[j] * dv[i];
        }
    }
}


}  // namespace reshala
