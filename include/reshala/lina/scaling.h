#pragma once

#include "reshala/lina/core/sparse_matrix.h"

namespace reshala {

struct Scaling {
    std::vector<Index> row;
    std::vector<Index> col;
};

struct ScaleReport {
    Scalar frob_norm = 0.0;
    Scalar max_min_ratio = 0.0;
};
std::ostream& operator<<(std::ostream& os, const ScaleReport& rep);

ScaleReport GetScaleReport(const SparseRowMatrix& A);

}  // namespace reshala
