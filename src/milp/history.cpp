#include "reshala/milp/history.h"

namespace reshala {

void History::Add(Direction dir, Index iv, Scalar dy, Scalar dx) {
    assert(!IsZero(MinFraction(dx)));
    assert(dx > 0);

    Scalar slope = std::max(dy, 0.0) / dx;  // Может быть околонулевой отрицательный мусор

    n_[iv][Dir2Index(dir)]++;
    s_[iv][Dir2Index(dir)] += slope;
    s2_[iv][Dir2Index(dir)] += slope * slope;
}

}  // namespace reshala
