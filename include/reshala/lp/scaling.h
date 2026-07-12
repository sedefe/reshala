#pragma once

#include "reshala/model/milp_model.h"

namespace reshala {

struct ScaleReport {
    ScaleReport() = default;
    ScaleReport(const SparseRowMatrix& A);

    Scalar frob_norm;
    Scalar max_min_ratio;
};
std::ostream& operator<<(std::ostream& os, const ScaleReport& rep);

struct ScaleStats {
    ScaleReport before;
    ScaleReport after;
};
std::ostream& operator<<(std::ostream& os, const ScaleStats& stats);

class Scaling {
   public:
    void ScaleModel(MilpModel& model);
    ScaleStats stats;

    std::vector<Index> row;
    std::vector<Index> col;
};

}  // namespace reshala
