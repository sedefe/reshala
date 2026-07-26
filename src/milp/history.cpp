#include "reshala/milp/history.h"

namespace reshala {

void History::Add(Index iv, Direction dir, Scalar dy, Scalar dx) {
    assert(!IsZero(MinFraction(dx)));
    assert(dx > 0);

    Scalar slope = std::max(dy, 0.0) / dx;  // Может быть околонулевой отрицательный мусор

    Index d = Dir2Index(dir);
    n_[iv][d]++;
    s_[iv][d] += slope;
    s2_[iv][d] += slope * slope;
}

Scalar History::Estimate(Index iv, Direction dir, Scalar dx) const {
    Index d = Dir2Index(dir);
    return s_[iv][d] / n_[iv][d];
}

Scalar History::GetSigma(Index iv, Direction dir) const {
    Index d = Dir2Index(dir);
    return (s2_[iv][d] - s_[iv][d] * s_[iv][d] / n_[iv][d]) / (n_[iv][d] - 1);
}

}  // namespace reshala
