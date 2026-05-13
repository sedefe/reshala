#include "reshala/linalg/sparse_matrix.h"

namespace reshala {

void Srm2Scm(const SparseRowMatrix& srm, SparseColMatrix& scm) {
    // Assuming scm is already resized
    auto m = srm.GetNRows();
    auto n = srm.GetNCols();
    for (Index i = 0; i < m; ++i) {
        const auto& row = srm.GetRows()[i];
        const auto& indices = row.indices();
        const auto& values = row.values();

        for (size_t j = 0; j < indices.size(); ++j) {
            auto col = indices[j];
            auto value = values[j];

            scm.GetCol(col).Push(i, value);
        }
    }
}

}  // namespace reshala
