#include "reshala/lina/core/operators.h"

#include <assert.h>

#include <cstring>

namespace reshala {

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

void dot(const SparseVector& sv1, const DenseVector& dv2, Scalar& res) {
    dot(dv2.data(), sv1, res);
}

void dot(const DenseVector& dv1, const DenseVector& dv2, Scalar& res) {
    assert(dv1.size() == dv2.size() && "Dv x Dv: incompatible sizes");

    res = 0;
    for (size_t i = 0; i < dv1.size(); ++i) {
        res += dv1[i] * dv2[i];
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

void MulDmSv(const DenseMatrix& dm, const SparseVector& sv, DenseVector& res) {
    auto m = dm.GetNRows();
    [[maybe_unused]] auto n = dm.GetNCols();
    assert(n == sv.dim() && "Dm x Sv: incompatible sizes");

    for (Index i = 0; i < m; i++) {
        const Scalar* row = dm[i];
        dot(row, sv, res[i]);
    }
}

bool Invert(DenseMatrix& dm) {
    const Scalar kPivotEps = 1e-12;
    const Index m = dm.GetNCols();
    assert(m == dm.GetNRows());

    // ---------- 1. In‑place LU decomposition ----------
    std::vector<Index> perm(m);
    for (Index i = 0; i < m; ++i) perm[i] = i;

    for (Index k = 0; k < m; ++k) {
        Scalar max_val = 0.0;
        Index pivot = k;
        for (Index i = k; i < m; ++i) {
            Scalar val = std::abs(dm[i][k]);
            if (val > max_val) {
                max_val = val;
                pivot = i;
            }
        }
        if (max_val < kPivotEps) return false;  // singular matrix

        if (pivot != k) {
            Scalar* row_k = dm[k];
            Scalar* row_p = dm[pivot];
            for (Index j = 0; j < m; ++j) std::swap(row_k[j], row_p[j]);
            std::swap(perm[k], perm[pivot]);
        }

        Scalar* row_k = dm[k];
        const Scalar diag = row_k[k];
        for (Index i = k + 1; i < m; ++i) {
            Scalar* row_i = dm[i];
            const Scalar mult = row_i[k] / diag;
            row_i[k] = mult;
            for (Index j = k + 1; j < m; ++j) {
                row_i[j] -= mult * row_k[j];
            }
        }
    }

    // ---------- 2. Solve for inverse columns, store in temporary buffer ----------
    std::vector<Scalar> inv(m * m);  // temporary storage for the inverse
    std::vector<Scalar> rhs(m), y(m), x(m);

    for (Index j = 0; j < m; ++j) {
        // Permuted identity column: P * e_j
        for (Index i = 0; i < m; ++i) {
            rhs[i] = (perm[i] == j) ? 1.0 : 0.0;
        }

        // Forward substitution: L y = rhs
        for (Index i = 0; i < m; ++i) {
            Scalar sum = 0.0;
            for (Index k = 0; k < i; ++k) {
                sum += dm[i][k] * y[k];
            }
            y[i] = rhs[i] - sum;
        }

        // Backward substitution: U x = y
        for (Index i = m - 1; i >= 0; --i) {
            Scalar sum = 0.0;
            for (Index k = i + 1; k < m; ++k) {
                sum += dm[i][k] * x[k];
            }
            const Scalar diag = dm[i][i];
            if (std::abs(diag) <= kPivotEps) return false;
            x[i] = (y[i] - sum) / diag;
        }

        // Store column j of the inverse in the temporary buffer
        for (Index i = 0; i < m; ++i) {
            inv[i * m + j] = x[i];
        }
    }

    // ---------- 3. Copy temporary buffer back into the original matrix ----------
    std::memcpy(dm[0], inv.data(), m * m * sizeof(Scalar));

    return true;
}

}  // namespace reshala
