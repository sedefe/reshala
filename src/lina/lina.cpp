#include "reshala/lina/lina.h"

#include "reshala/lina/core/operators.h"

namespace reshala {

Lina::Lina(const SparseColMatrix& Ac, const SparseRowMatrix& Ar, const LpBasis* basis)
    : Ac_(Ac),
      Ar_(Ar),
      m(Ac.GetNRows()),
      n(Ac.GetNCols()),
      basis_(basis),
      Lr(m, m),
      Ur(m, m),
      Lc(m, m),
      Uc(m, m),
      u_diag(m),
      ftran_res(m) {
    Ar_.Resize(m, n + m);
    Ac_.Resize(m, n + m);
    for (Index i = 0; i < m; i++) {
        Ar_.GetRow(i).Push(n + i, 1.0);
        Ac_.GetCol(n + i) = {m, i, 1.0};
    }

    Scale();
    Reset();
}

std::ostream& operator<<(std::ostream& os, const LinaStats& stats) {
    os << "Lina: " << stats.n_lus << " LUs, " << stats.n_updates << " updates\n";
    return os;
}

void Lina::Reset() {
    Binv_.ResizeAsZero(m, m);
    for (Index i = 0; i < m; i++) {
        Binv_[i][i] = 1;
    }

    row_perm.resize(m);
    row_perm_inv.resize(m);
    Lr.Clear();
    Ur.Clear();
    Lc.Clear();
    Uc.Clear();
    for (Index i = 0; i < m; ++i) {
        row_perm[i] = row_perm_inv[i] = i;
        u_diag[i] = std::ldexp(1.0, -scaling.row[i]);
    };
}

void Lina::Scale() {
    scaling.row.assign(m, 0.0);
    scaling.col.assign(n + m, 0.0);

    // ScaleReport rep0 = GetScaleReport(Ar_);
    // std::cout << rep0;

    // Init scaling
    Index exp;
    Index sum;
    Index kMaxIters = 5;

    for (Index i = 0; i < kMaxIters; i++) {
        for (Index ic = 0; ic < m; ic++) {
            if (Ar_.GetRow(ic).Size() == 0) continue;
            sum = 0;
            for (SvIterator el(Ar_.GetRow(ic)); el; ++el) {
                std::frexp(el.value(), &exp);
                sum += exp - scaling.col[el.index()];
            }
            scaling.row[ic] = Scalar(sum) / Ar_.GetRow(ic).Size() + 0.5;
        }

        for (Index iv = 0; iv < n; iv++) {
            if (Ac_.GetCol(iv).Size() == 0) continue;
            sum = 0;
            for (SvIterator el(Ac_.GetCol(iv)); el; ++el) {
                std::frexp(el.value(), &exp);
                sum += exp - scaling.row[el.index()];
            }
            scaling.col[iv] = Scalar(sum) / Ac_.GetCol(iv).Size() + 0.5;
        }
    }

    // Apply
    for (Index iv = 0; iv < n + m; iv++) {
        for (MutableSvIterator el(Ac_.GetCol(iv)); el; ++el) {
            el.valueRef() = std::ldexp(el.value(), -scaling.row[el.index()] - scaling.col[iv]);
        }
    }
    for (Index ic = 0; ic < m; ic++) {
        for (MutableSvIterator el(Ar_.GetRow(ic)); el; ++el) {
            el.valueRef() = std::ldexp(el.value(), -scaling.row[ic] - scaling.col[el.index()]);
        }
    }

    // ScaleReport rep1 = GetScaleReport(Ar_);
    // std::cout << rep1;
}

}  // namespace reshala
