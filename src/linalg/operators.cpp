#include "reshala/linalg/operators.h"

#include <assert.h>

namespace reshala {

SparseVector sub(const SparseVector& sv1, const SparseVector& sv2) {
    assert(sv1.dim() == sv2.dim() && "SparseVector sub: vectors are of different dimentions");

    std::vector<Index> result_indices;
    std::vector<Scalar> result_values;

    auto& x_indices = sv1.indices();
    auto& x_values = sv1.values();
    auto& y_indices = sv2.indices();
    auto& y_values = sv2.values();

    size_t i = 0, j = 0;
    size_t x_size = x_indices.size();
    size_t y_size = y_indices.size();

    while (i < x_size || j < y_size) {
        Index idx;
        Scalar x_val = Scalar(0);
        Scalar y_val = Scalar(0);

        if (i < x_size && (j >= y_size || x_indices[i] < y_indices[j])) {
            idx = x_indices[i];
            x_val = x_values[i];
            i++;
        } else if (j < y_size && (i >= x_size || y_indices[j] < x_indices[i])) {
            idx = y_indices[j];
            y_val = y_values[j];
            j++;
        } else {
            idx = x_indices[i];
            x_val = x_values[i];
            y_val = y_values[j];
            i++;
            j++;
        }

        Scalar result_val = x_val - y_val;

        if (!IsZero(result_val)) {
            result_indices.push_back(idx);
            result_values.push_back(result_val);
        }
    }

    return SparseVector(sv1.dim(), result_indices, result_values);
}

void dot(const SparseVector& sv1, const SparseVector& sv2, Scalar& res) {
    assert(sv1.dim() == sv2.dim() && "SparseVector dot: vectors are of different dimensions");

    auto& x_indices = sv1.indices();
    auto& x_values = sv1.values();
    auto& y_indices = sv2.indices();
    auto& y_values = sv2.values();

    Index i = 0, j = 0;
    auto size1 = sv1.Size();
    auto size2 = sv2.Size();

    res = 0;

    while (i < size1 && j < size2) {
        if (x_indices[i] < y_indices[j]) {
            i++;
        } else if (y_indices[j] < x_indices[i]) {
            j++;
        } else {
            res += x_values[i] * y_values[j];
            i++;
            j++;
        }
    }
}

void dot(const Scalar* ar1, const SparseVector& sv2, Scalar& res) {
    res = 0;
    const auto& indices = sv2.indices();
    const auto& values = sv2.values();

    for (size_t i = 0; i < indices.size(); ++i) {
        res += ar1[indices[i]] * values[i];
    }
}

void dot(const DenseVector& dv1, const SparseVector& sv2, Scalar& res) {
    assert(dv1.size() == sv2.dim() && "Dv x Sv: incompatible sizes");
    dot(dv1.data(), sv2, res);
}

void dot(const DenseVector& dv1, const DenseVector& dv2, Scalar& res) {
    assert(dv1.size() == dv2.size() && "Dv x Dv: incompatible sizes");

    res = 0;
    for (size_t i = 0; i < dv1.size(); ++i) {
        res += dv1[i] * dv2[i];
    }
}

void MulScmSv(const SparseColMatrix& scm, const SparseVector& sv, DenseVector& res) {
    auto m = scm.GetNRows();
    auto n = scm.GetNCols();
    assert(n == sv.dim() && "Scm x Sv: incompatible sizes");

    res.assign(m, Scalar(0));
    for (SvIterator el_i(sv); el_i; ++el_i) {
        const auto& col = scm.GetCol(el_i.index());
        for (SvIterator el_j(col); el_j; ++el_j) {
            res[el_j.index()] += el_i.value() * el_j.value();
        }
    }
}

void MulScmDv(const SparseColMatrix& scm, const DenseVector& dv, DenseVector& res) {
    auto m = scm.GetNRows();
    auto n = scm.GetNCols();
    assert(n == dv.size() && "Scm x Dv: incompatible sizes");

    res.assign(m, Scalar(0));
    for (Index i = 0; i < n; i++) {
        if (IsZero(dv[i])) continue;
        const auto& col = scm.GetCol(i);
        for (SvIterator el(col); el; ++el) {
            res[el.index()] += dv[i] * el.value();
        }
    }
}

void MulDvDm(const DenseVector& dv, const DenseMatrix& dm, DenseVector& res) {
    auto m = dm.GetNRows();
    auto n = dm.GetNCols();
    assert(m == dv.size() && "Dv x Dm: incompatible sizes");

    res.assign(n, 0.0);
    for (Index i = 0; i < m; i++) {
        const Scalar* row = dm.RowView(i);
        for (Index j = 0; j < n; j++) {
            res[j] += row[j] * dv[i];
        }
    }
}

void MulDmDv(const DenseMatrix& dm, const DenseVector& dv, DenseVector& res) {
    auto m = dm.GetNRows();
    auto n = dm.GetNCols();
    assert(n == dv.size() && "Dm x Dv: incompatible sizes");

    res.assign(m, 0.0);
    for (Index i = 0; i < m; i++) {
        const Scalar* row = dm.RowView(i);
        Scalar sum = 0.0;
        for (Index j = 0; j < n; j++) {
            sum += row[j] * dv[j];
        }
        res[i] = sum;
    }
}

void MulDmSv(const DenseMatrix& dm, const SparseVector& sv, DenseVector& res) {
    auto m = dm.GetNRows();
    auto n = dm.GetNCols();
    assert(n == sv.dim() && "Dm x Sv: incompatible sizes");

    for (Index i = 0; i < m; i++) {
        const Scalar* row = dm.RowView(i);
        dot(row, sv, res[i]);
    }
}

}  // namespace reshala
