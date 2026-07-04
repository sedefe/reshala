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
      ftran_res(m) {
    Ar_.Resize(m, n + m);
    Ac_.Resize(m, n + m);
    for (Index i = 0; i < m; i++) {
        Ar_.GetRow(i).Push(n + i, 1.0);
        Ac_.GetCol(n + i) = {m, i, 1.0};
    }
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
    Lc.Clear();
    for (Index i = 0; i < m; ++i) {
        row_perm[i] = row_perm_inv[i] = i;
        Ur.GetRow(i) = SparseVector(m, i, 1.0);
        Uc.GetCol(i) = SparseVector(m, i, 1.0);
    };
}

}  // namespace reshala
