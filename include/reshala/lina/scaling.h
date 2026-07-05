#pragma once

#include "reshala/lina/core/sparse_matrix.h"

namespace reshala {

struct Scaling {
    std::vector<Index> row;
    std::vector<Index> col;
};

struct ScaleReport {
    Scalar frob_norm;
    Scalar max_min_ratio;
};
std::ostream& operator<<(std::ostream& os, const ScaleReport& rep);

ScaleReport GetScaleReport(const SparseRowMatrix& A);

}  // namespace reshala
