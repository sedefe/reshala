#include "reshala/lina/lina.h"

#include "reshala/lina/core/operators.h"

namespace reshala {

void Lina::Reset() {
    Binv_.ResizeAsZero(m, m);
    for (Index iv = 0; iv < m; iv++) {
        Binv_[iv][iv] = 1;
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
