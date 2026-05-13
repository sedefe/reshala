#include "reshala/presolve/model_info.h"

namespace reshala {

void ModelInfo::CompressVars() {
    auto m = model_.GetNCons();
    auto n = model_.GetNVars();

    // Objective & Domain
    auto& coeffs = model_.GetObj().coefficients;
    auto& domain = model_.GetDomain();

    Index i_write = 0;
    for (Index i_read = 0; i_read < n; ++i_read) {
        if (!var_mask.Get(i_read)) {
            if (i_write != i_read) {
                coeffs[i_write] = std::move(coeffs[i_read]);
                model_.GetAc().GetCol(i_write) = std::move(model_.GetAc().GetCol(i_read));
                domain.Move(i_read, i_write);
            }
            ++i_write;
        }
    }

    // Mapping from old column to new column index
    std::vector<Index> new_index_map(n, -1);
    Index next_new_index = 0;

    for (Index iv = 0; iv < n; ++iv) {
        if (!var_mask.Get(iv)) {
            new_index_map[iv] = next_new_index++;
        }
    }

    for (SparseVector& row : model_.GetAr().GetRows()) {
        Index i_write = 0;

        for (Index i_read = 0; i_read < row.indices().size(); ++i_read) {
            Index new_col = new_index_map[row.indices()[i_read]];

            if (new_col != -1) {
                if (i_write != i_read) {
                    row.indices()[i_write] = new_col;
                    row.values()[i_write] = row.values()[i_read];
                } else {
                    row.indices()[i_write] = new_col;
                }
                ++i_write;
            }
        }

        row.Resize(i_write);
    }

    model_.Resize(m, n - deleted_vars.size());

    con_mask.Clear();
    var_mask.Clear();
    deleted_cons.clear();
    deleted_vars.clear();
}

}  // namespace reshala
