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

    Refactor();
}

std::ostream& operator<<(std::ostream& os, const LinaStats& stats) {
    Scalar avg_fillin = stats.total_nnz_b > 0
                            ? Scalar(stats.total_nnz_l + stats.total_nnz_u) / stats.total_nnz_b
                            : 0.0;
    os << "Lina stats: \n"
       << "\t" << stats.n_lus << " LUs, " << stats.n_updates << " updates\n"
       << "\tAvg fill-in: " << avg_fillin << "\n";
    return os;
}

}  // namespace reshala
