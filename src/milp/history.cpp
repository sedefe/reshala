#include "reshala/milp/history.h"

namespace reshala {

void History::Add(Index iv, Direction dir, Scalar dy, Scalar dx) {
    assert(!IsZero(MinFraction(dx)));
    assert(StrongGt(dx, 0.0));

    Scalar slope = std::max(dy, 0.0) / dx;  // Может быть околонулевой отрицательный мусор

    Index d = Dir2Index(dir);
    n_[iv][d]++;
    s_[iv][d] += slope;
    s2_[iv][d] += slope * slope;
}

bool History::IsEnough(Index iv) const {
    // m >= 2*s  =>  m*m >= 4*v

    if (n_[iv][0] < kMinSamples or n_[iv][1] < kMinSamples) return false;

    Scalar m0 = GetMean(iv, Direction::kLeft);
    Scalar m1 = GetMean(iv, Direction::kRight);
    Scalar v0 = GetVariance(iv, Direction::kLeft);
    Scalar v1 = GetVariance(iv, Direction::kRight);

    return (m0 * m0 >= 4 * v0) and (m1 * m1 >= 4 * v1);
}

}  // namespace reshala
