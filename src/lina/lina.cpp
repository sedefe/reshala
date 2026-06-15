#include "reshala/lina/lina.h"

#include "reshala/lina/core/operators.h"

namespace reshala {

void Lina::Init() {
    Binv_.ResizeAsZero(m, m);
    for (Index iv = 0; iv < m; iv++) {
        Binv_[iv][iv] = 1;
    }
    n_updates_ = 0;

    perm.resize(m);
    Lr.Clear();
    Ur.Clear();
    Lc.Clear();
    Uc.Clear();
    for (Index i = 0; i < m; ++i) {
        perm[i] = i;
        Ur.GetRow(i).Push(i, 1.0);
        Uc.GetCol(i).Push(i, 1.0);
    };
}

}  // namespace reshala
