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
    std::vector<Index> new_index_map(n, -1);
    // Todo: use (sorted) deleted_vars_
    for (Index i_read = 0; i_read < n; ++i_read) {
        if (!GetVarMask(i_read)) {
            if (i_write != i_read) {
                orig_var_idx_[i_write] = std::move(orig_var_idx_[i_read]);
                coeffs[i_write] = std::move(coeffs[i_read]);
                model_.GetCol(i_write) = std::move(model_.GetCol(i_read));
                domain.Move(i_read, i_write);
            }
            new_index_map[i_read] = i_write;
            ++i_write;
        }
    }

    // Ar
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

    // Implications
    Index new_n_impl = 0;
    for (Index i = 0; i < implications_.size(); i++) {
        Index iv1 = implications_[i].x_ind;
        Index iv2 = implications_[i].y_ind;
        if (GetVarMask(iv1) or GetVarMask(iv2)) continue;
        implications_[i].x_ind = new_index_map[iv1];
        implications_[i].y_ind = new_index_map[iv2];
        if (new_n_impl != i) std::swap(implications_[i], implications_[new_n_impl]);
        new_n_impl++;
    }
    implications_.resize(new_n_impl);

    model_.Resize(m, n - deleted_vars_.size());

    var_mask_.Clear();
    deleted_vars_.clear();
}

void ModelTracker::CalcActivities() {
    auto m = model_.GetNCons();
    activities_.resize(m);
    for (Index ic = 0; ic < m; ic++) {
        activities_[ic] = CalcActivity(ic);
    }
}

Activity ModelTracker::CalcActivity(Index ic) const {
    Activity act;
    for (SvIterator el(model_.GetRow(ic)); el; ++el) {
        act.AddTerm(el.value(), model_.GetBounds(el.index()));
    }
    return act;
}

void ModelTracker::FixVar(Index iv, Scalar val) {
    model_.GetObj().c0 += model_.GetObj().coefficients[iv] * val;

    UpdVarBounds(iv, {0, 0});  // Убираем эту переменную из активити

    for (SvIterator el(model_.GetCol(iv)); el; ++el) {
        const Bounds& rhs = model_.GetRhs(el.index());
        model_.GetRhs(el.index()) = {rhs.le - el.value() * val, rhs.ri - el.value() * val};
    }

    transforms_.push_back(
        std::make_unique<FixVariableTransform>(FixVariableTransform(orig_var_idx_[iv], val)));
    MaskVar(iv);
}

bool ModelTracker::SimpleSub(Index iv1, Scalar a, Index iv2, Scalar b) {
    // iv1 <- a*iv2 + b
    const Bounds& bnd1 = model_.GetBounds(iv1);
    const Bounds& bnd2 = model_.GetBounds(iv2);

    // Границы iv1 могут повлиять на границы iv2
    Bounds derived_bnd2 = (a >= 0) ? Bounds((bnd1.le - b) / a, (bnd1.ri - b) / a)
                                   : Bounds((bnd1.ri - b) / a, (bnd1.le - b) / a);
    if (model_.GetIntegrality(iv2)) {
        derived_bnd2.le = WeakCeil(derived_bnd2.le);
        derived_bnd2.ri = WeakFloor(derived_bnd2.ri);
    }

    const Bounds new_bnd2 = BoundsIntersection(bnd2, derived_bnd2);
    if (StrongGt(new_bnd2.le, new_bnd2.ri)) return false;

    model_.GetObj().c0 += model_.GetObj().coefficients[iv1] * b;
    model_.GetObj().coefficients[iv2] += a * model_.GetObj().coefficients[iv1];

    for (SvIterator el(model_.GetCol(iv1)); el; ++el) {
        Index ic = el.index();
        if (GetConMask(ic)) continue;

        Scalar val_iv1 = el.value();
        Scalar val_iv2 = a * val_iv1;
        SparseVector& row = model_.GetRow(el.index());
        row.EraseIndex(iv1);

        Activity& act = activities_[ic];
        act.RmTerm(val_iv1, bnd1);

        // Ar & activity
        Index offset = row.FindOffset(iv2);  // Todo double iteration on iv1&iv2?
        if (row.indices()[offset] == iv2) {  // Can update inplace
            // Scalar
            Scalar old_val_iv2 = row.values()[offset];
            Scalar new_val_iv2 = old_val_iv2 + val_iv2;
            row.values()[offset] = new_val_iv2;
            // Проверяем на 0, иначе ниже при апдейте столбца Ac и Ar станут неконсистентны
            if (IsZero(row.values()[offset])) {
                row.EraseOffset(offset);
            }
            act.RmTerm(old_val_iv2, bnd2);
            act.AddTerm(new_val_iv2, bnd2);
        } else {  // Insert a new value
            if (!IsZero(val_iv2)) {
                row.Insert(iv2, val_iv2, offset);
            }
            act.AddTerm(val_iv2, bnd2);
        }

        // Rhs
        const Bounds& rhs = model_.GetRhs(el.index());
        model_.GetRhs(el.index()) = {rhs.le - val_iv1 * b, rhs.ri - val_iv1 * b};
    }
    // Ac
    model_.GetCol(iv2) = model_.GetCol(iv2) + a * model_.GetCol(iv1);

    if (StrongGt(new_bnd2.le, bnd2.le) or StrongLt(new_bnd2.ri, bnd2.ri)) {
        UpdVarBounds(iv2, new_bnd2);
    }

    transforms_.push_back(std::make_unique<SimpleSubTransform>(
        SimpleSubTransform(orig_var_idx_[iv1], a, orig_var_idx_[iv2], b)));
    MaskVar(iv1);

    return true;
}

void ModelTracker::SlackSub(Index ic, Index iv, Scalar a) {
    // Sum(a_k*iv_k) + a * iv = b
    const Bounds bnd = model_.GetBounds(iv);
    Bounds& rhs = model_.GetRhs(ic);
    Scalar b = (rhs.le + rhs.ri) / 2;

    if (a >= 0) {
        rhs.le = b - a * bnd.ri;
        rhs.ri = b - a * bnd.le;
    } else {
        rhs.le = b - a * bnd.le;
        rhs.ri = b - a * bnd.ri;
    }

    Scalar c = model_.GetObj().coefficients[iv];
    SparseVector sv(model_.GetNVars());
    model_.GetObj().c0 += model_.GetObj().coefficients[iv] * b / a;
    for (SvIterator el(model_.GetRow(ic)); el; ++el) {
        if (el.index() == iv) continue;
        if (GetVarMask(el.index())) continue;
        model_.GetObj().coefficients[el.index()] -= c * el.value() / a;
        sv.Push(orig_var_idx_[el.index()], -el.value() / a);
    }

    UpdVarBounds(iv, {0., 0.});

    transforms_.push_back(
        std::make_unique<LinCombTransform>(LinCombTransform(orig_var_idx_[iv], sv, b / a)));
    MaskVar(iv);
}

void ModelTracker::UpdRhs(Index ic, Bounds rhs) {
    if (IsZero(rhs.le - rhs.ri)) rhs.le = rhs.ri = (rhs.le + rhs.ri) / 2;

    model_.GetRhs(ic) = rhs;
    stat.n_ch_rhs++;
}

void ModelTracker::UpdVarBounds(Index iv, Bounds bnd) {
    const Bounds& old_bnd = model_.GetBounds(iv);
    if (IsZero(bnd.le - bnd.ri)) bnd.le = bnd.ri = (bnd.le + bnd.ri) / 2;

    for (SvIterator el(model_.GetCol(iv)); el; ++el) {
        activities_[el.index()].RmTerm(el.value(), old_bnd);
        activities_[el.index()].AddTerm(el.value(), bnd);
    }
    model_.SetBounds(iv, bnd);
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
    activities_[ic].RmTerm(value_ref, bnd);
    activities_[ic].AddTerm(val, bnd);

    // Coeffs
    model_.GetCol(iv).AtRef(ic) = val;
    value_ref = val;
    stat.n_ch_coeff++;
}

void ModelTracker::ScaleObj(Scalar x) {
    for (auto& val : model_.GetObj().coefficients) val *= x;
    model_.GetObj().c0 *= x;
    model_.GetObj().mult /= x;
}

void ModelTracker::ScaleObjExp(Index e) {
    for (auto& val : model_.GetObj().coefficients) val = std::ldexp(val, e);
    model_.GetObj().c0 = std::ldexp(model_.GetObj().c0, e);
    model_.GetObj().mult = std::ldexp(model_.GetObj().mult, -e);
}

void ModelTracker::ScaleRow(Index ic, Scalar x) {
    SparseVector& row = model_.GetRow(ic);
    row *= x;
    for (SvIterator el(row); el; ++el) {
        auto iv = el.index();
        if (GetVarMask(iv)) continue;
        model_.GetCol(iv).AtRef(ic) *= x;
    }
    activities_[ic].Scale(x);

    const auto& rhs = model_.GetRhs(ic);
    UpdRhs(ic, {rhs.le * x, rhs.ri * x});

    stat.n_ch_coeff++;
}

Bounds ModelTracker::DeriveBounds(Index ic, Index iv, Activity act, const Bounds& bnd,
                                  Scalar val) const {
    // Todo: pass non-const bounds and modify them
    // We will modify the activity, so it is passed as a value
    Bounds derived;
    const Bounds& rhs = model_.GetRhs(ic);
    act.RmTerm(val, bnd);
    auto lhs = act.GetRange();
    if (val > 0) {
        derived.le = (rhs.le - lhs.ri) / val;
        derived.ri = (rhs.ri - lhs.le) / val;
    } else {
        derived.le = (rhs.ri - lhs.le) / val;
        derived.ri = (rhs.le - lhs.ri) / val;
    }
    if (model_.GetIntegrality(iv)) {
        derived.le = WeakCeil(derived.le);
        derived.ri = WeakFloor(derived.ri);
    }
    return derived;
}

}  // namespace reshala
