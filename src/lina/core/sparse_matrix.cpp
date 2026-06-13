#include "reshala/lina/core/sparse_matrix.h"

namespace reshala {

void Srm2Scm(const SparseRowMatrix& srm, SparseColMatrix& scm) {
    // Assuming scm is already resized
    scm.Clear();
    auto m = srm.GetNRows();
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

std::ostream& operator<<(std::ostream& os, const SparseRowMatrix& svm) {
    Index m = svm.GetNRows();
    for (Index i = 0; i < m; i++) {
        os << "row " << i << ": " << svm.GetRow(i) << "\n";
    }
    return os;
}

}  // namespace reshala
