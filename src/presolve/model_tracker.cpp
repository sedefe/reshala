#include "reshala/presolve/model_tracker.h"

namespace reshala {

ModelTracker::ModelTracker(MilpModel& model)
    : model_(model), con_mask_(model.GetNCons()), var_mask_(model.GetNVars()) {
    orig_n_vars_ = model.GetNVars();
    orig_var_idx_.resize(orig_n_vars_);
    for (Index iv = 0; iv < orig_n_vars_; ++iv) {
        orig_var_idx_[iv] = iv;
    }
}

void ModelTracker::CompressCons() {
    auto m = model_.GetNCons();
    auto n = model_.GetNVars();

    // Ar & Rhs & Activities
    auto& rhs = model_.GetRhs();

    Index i_write = 0;
    // Todo: use (sorted) deleted_cons_
    for (Index i_read = 0; i_read < m; ++i_read) {
        if (!con_mask_.Get(i_read)) {
            if (i_write != i_read) {
                model_.GetRow(i_write) = std::move(model_.GetRow(i_read));
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

void ModelTracker::CompressVars() {
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
                model_.GetCol(i_write) = std::move(model_.GetCol(i_read));
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

void ModelTracker::CalcActivities() {
    auto m = model_.GetNCons();
    activities_.resize(m);
    for (Index ic = 0; ic < m; ic++) {
        Bounds act = {0, 0};
        for (SvIterator el(model_.GetRow(ic)); el; ++el) {
            const Bounds& bnd = model_.GetBounds(el.index());
            Scalar val = el.value();
            if (val >= 0) {
                act.le += val * bnd.le;
                act.ri += val * bnd.ri;
            } else {
                act.le += val * bnd.ri;
                act.ri += val * bnd.le;
            }
        }
        activities_[ic] = act;
    }
}

void ModelTracker::FixVar(Index iv, Scalar val) {
    model_.GetObj().c0 += model_.GetObj().mult * (model_.GetObj().coefficients[iv] * val);

    UpdVarBounds(iv, {0, 0});  // Убираем эту переменную из активити

    for (SvIterator el(model_.GetCol(iv)); el; ++el) {
        const Bounds& rhs = model_.GetRhs(el.index());
        model_.GetRhs(el.index()) = {rhs.le - el.value() * val, rhs.ri - el.value() * val};
    }

    transforms_.push_back(
        std::make_unique<FixVariableTransform>(FixVariableTransform(orig_var_idx_[iv], val)));
    MaskVar(iv);
}

void ModelTracker::SimpleSub(Index iv1, Scalar a, Index iv2, Scalar b) {
    // iv1 <- a*iv2 + b
    model_.GetObj().c0 += model_.GetObj().mult * (model_.GetObj().coefficients[iv1] * b);
    model_.GetObj().coefficients[iv2] += a * model_.GetObj().coefficients[iv1];

    for (SvIterator el(model_.GetCol(iv1)); el; ++el) {
        if (GetConMask(el.index())) continue;

        // Todo: double iteration on col[iv1] & col[iv2]
        Scalar val = el.value();
        SparseVector& row = model_.GetRow(el.index());
        row.Erase(iv1);

        auto pos = row.FindIndex(iv2);
        if (*pos == iv2) {  // Can update inplace
            row.values()[pos - row.indices().begin()] += a * val;
            // Проверяем на 0, иначе ниже при апдейте столбца Ac и Ar станут неконсистентны
            if (IsZero(row.values()[pos - row.indices().begin()])) {
                row.Erase(pos);
            }
        } else {  // Insert a new value
            if (!IsZero(a * val)) {
                row.Insert(iv2, a * val, pos);
            }
        }

        const Bounds& rhs = model_.GetRhs(el.index());
        model_.GetRhs(el.index()) = {rhs.le - el.value() * b, rhs.ri - el.value() * b};
    }
    model_.GetCol(iv2) = model_.GetCol(iv2) + a * model_.GetCol(iv1);

    transforms_.push_back(std::make_unique<SimpleSubTransform>(
        SimpleSubTransform(orig_var_idx_[iv1], a, orig_var_idx_[iv2], b)));
    MaskVar(iv1);
}

void ModelTracker::UpdRhs(Index ic, const Bounds& bnd) {
    model_.GetRhs(ic) = bnd;
    stat.n_ch_rhs++;
}

void ModelTracker::UpdVarBounds(Index iv, const Bounds& bnd) {
    const Bounds& old_bnd = model_.GetBounds(iv);
    const Bounds diff = {bnd.le - old_bnd.le, bnd.ri - old_bnd.ri};

    for (SvIterator el(model_.GetCol(iv)); el; ++el) {
        const Bounds& act = activities_[el.index()];
        Scalar val = el.value();
        activities_[el.index()] =
            (val >= 0)  // Todo Логика повторяется в 3.2 и 3.3. Перенести в утилиты?
                ? Bounds{act.le + val * diff.le, act.ri + val * diff.ri}
                : Bounds{act.le + val * diff.ri, act.ri + val * diff.le};
    }
    model_.SetBounds(iv, bnd);

    if (model_.GetType(iv) == BndType::kInfeasible) {
        infeasible_ = true;
    }
    stat.n_ch_bnd++;
}

void ModelTracker::UpdCoeff(Index ic, Index iv, Scalar val) {
    // Todo Вообще это можно быстрее делать, т.к. к части из следующих операций доступ есть из тех
    // мест, где мы вызываем этот апдейт
    // Todo handle near-zeros
    Scalar& value_ref = model_.GetRow(ic).AtRef(iv);
    assert((value_ref >= 0 and val >= 0) or (value_ref <= 0 and val <= 0));
    auto d = val - value_ref;
    const Bounds& bnd = model_.GetBounds(iv);

    if (val >= 0) {
        activities_[ic].ri -= d * bnd.ri;
        activities_[ic].le -= d * bnd.le;
    } else {
        activities_[ic].ri -= d * bnd.le;
        activities_[ic].le -= d * bnd.ri;
    }

    // Coeffs
    model_.GetCol(iv).AtRef(ic) = val;
    value_ref = val;
    stat.n_ch_coeff++;
}

void ModelTracker::ScaleObj(Scalar x) {
    for (auto& val : model_.GetObj().coefficients) val *= x;
    model_.GetObj().mult /= x;
}

void ModelTracker::ScaleRow(Index ic, Scalar x) {
    SparseVector& row = model_.GetRow(ic);
    row *= x;
    for (SvIterator el(row); el; ++el) {
        auto iv = el.index();
        if (GetVarMask(iv)) continue;
        model_.GetCol(iv).AtRef(ic) *= x;
    }
    activities_[ic].le *= x;
    activities_[ic].ri *= x;

    const auto& rhs = model_.GetRhs(ic);
    UpdRhs(ic, {rhs.le * x, rhs.ri * x});

    stat.n_ch_coeff++;
}

}  // namespace reshala
