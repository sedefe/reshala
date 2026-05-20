#include "reshala/presolve/model_info.h"

namespace reshala {

ModelInfo::ModelInfo(MilpModel& model)
    : model_(model), con_mask_(model.GetNCons()), var_mask_(model.GetNVars()) {
    orig_n_vars_ = model.GetNVars();
    orig_var_idx_.resize(orig_n_vars_);
    for (Index iv = 0; iv < orig_n_vars_; ++iv) {
        orig_var_idx_[iv] = iv;
    }
}

void ModelInfo::CompressCons() {
    auto m = model_.GetNCons();
    auto n = model_.GetNVars();

    // Ar & Rhs & Activities
    auto& rhs = model_.GetRhs();

    Index i_write = 0;
    // Todo: use (sorted) deleted_cons_
    for (Index i_read = 0; i_read < m; ++i_read) {
        if (!con_mask_.Get(i_read)) {
            if (i_write != i_read) {
                model_.GetAr().GetRow(i_write) = std::move(model_.GetAr().GetRow(i_read));
                rhs[i_write] = std::move(rhs[i_read]);
                activities_[i_write] = std::move(activities_[i_read]);
            }
            ++i_write;
        }
    }

    // Ac
    std::vector<Index> new_index_map(m, -1);
    Index next_new_index = 0;
    // Todo: use (sorted) deleted_cons_
    for (Index ic = 0; ic < m; ++ic) {
        if (!con_mask_.Get(ic)) {
            new_index_map[ic] = next_new_index++;
        }
    }

    for (SparseVector& col : model_.GetAc().GetCols()) {
        Index i_write = 0;
        for (Index i_read = 0; i_read < col.indices().size(); ++i_read) {
            Index new_col = new_index_map[col.indices()[i_read]];
            if (new_col != -1) {
                if (i_write != i_read) {
                    col.indices()[i_write] = new_col;
                    col.values()[i_write] = col.values()[i_read];
                } else {
                    col.indices()[i_write] = new_col;
                }
                ++i_write;
            }
        }

        col.Resize(i_write);
    }

    model_.Resize(m - deleted_cons_.size(), n);
    activities_.resize(m - deleted_cons_.size());

    con_mask_.Clear();
    deleted_cons_.clear();
}

void ModelInfo::CompressVars() {
    auto m = model_.GetNCons();
    auto n = model_.GetNVars();

    // Objective & Domain & Ac
    auto& coeffs = model_.GetObj().coefficients;
    auto& domain = model_.GetDomain();

    Index i_write = 0;
    // Todo: use (sorted) deleted_vars_
    for (Index i_read = 0; i_read < n; ++i_read) {
        if (!var_mask_.Get(i_read)) {
            if (i_write != i_read) {
                orig_var_idx_[i_write] = std::move(orig_var_idx_[i_read]);
                coeffs[i_write] = std::move(coeffs[i_read]);
                model_.GetAc().GetCol(i_write) = std::move(model_.GetAc().GetCol(i_read));
                domain.Move(i_read, i_write);
            }
            ++i_write;
        }
    }

    // Ar
    std::vector<Index> new_index_map(n, -1);
    Index next_new_index = 0;
    // Todo: use (sorted) deleted_vars_
    for (Index iv = 0; iv < n; ++iv) {
        if (!var_mask_.Get(iv)) {
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

    model_.Resize(m, n - deleted_vars_.size());

    var_mask_.Clear();
    deleted_vars_.clear();
}

void ModelInfo::CalcActivities() {
    auto m = model_.GetNCons();
    activities_.resize(m);
    for (Index ic = 0; ic < m; ic++) {
        Bounds act = {0, 0};
        for (SvIterator el(model_.GetAr().GetRow(ic)); el; ++el) {
            const Bounds& bnd = model_.GetBounds(el.index());
            if (el.value() >= 0) {
                act.le += el.value() * bnd.le;
                act.ri += el.value() * bnd.ri;
            } else {
                act.le += el.value() * bnd.ri;
                act.ri += el.value() * bnd.le;
            }
        }
        activities_[ic] = act;
    }
}

void ModelInfo::FixVar(Index iv, Scalar val) {
    model_.GetObj().c0 += model_.GetObj().mult * (model_.GetObj().coefficients[iv] * val);

    UpdVarBounds(iv, {0, 0});  // Убираем эту переменную из активити

    for (SvIterator el(model_.GetAc().GetCol(iv)); el; ++el) {
        const Bounds& rhs = model_.GetRhs(el.index());
        model_.GetRhs(el.index()) = {rhs.le - el.value() * val, rhs.ri - el.value() * val};
    }
}

void ModelInfo::UpdVarBounds(Index iv, const Bounds& bnd) {
    const Bounds& old_bnd = model_.GetBounds(iv);
    const Bounds diff = {bnd.le - old_bnd.le, bnd.ri - old_bnd.ri};

    for (SvIterator el(model_.GetAc().GetCol(iv)); el; ++el) {
        const Bounds& act = activities_[el.index()];
        activities_[el.index()] =
            (el.value() >= 0)
                ? Bounds{act.le + el.value() * diff.le, act.ri + el.value() * diff.ri}
                : Bounds{act.le + el.value() * diff.ri, act.ri + el.value() * diff.le};
    }
    model_.SetBounds(iv, bnd);

    if (model_.GetType(iv) == BndType::kInfeasible) {
        infeasible_ = true;
    }
}

}  // namespace reshala
