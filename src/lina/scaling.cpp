#include "reshala/lina/scaling.h"

namespace reshala {

std::ostream& operator<<(std::ostream& os, const ScaleReport& rep) {
    os << "FrobNorm: " << rep.frob_norm << ", MaxMinRatio: " << rep.max_min_ratio << "\n";
    return os;
}

ScaleReport GetScaleReport(const SparseRowMatrix& A) {
    ScaleReport rep;
    rep.frob_norm = 0.0;
    Index nnz = A.GetNnz();
    Scalar mx = -kInf, mn = kInf, val;
    for (Index ic = 0; ic < A.GetNRows(); ic++) {
        if (A.GetRow(ic).Size() == 0) continue;
        for (SvIterator el(A.GetRow(ic)); el; ++el) {
            val = std::abs(el.value());
            rep.frob_norm += val * val;
            mx = std::max(mx, val);
            mn = std::min(mn, val);
        }
    }
    rep.max_min_ratio = mx / mn;
    return rep;
}

}  // namespace reshala
