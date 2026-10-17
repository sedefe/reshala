#include "reshala/presolve/model_tracker.h"

namespace reshala {

ModelTracker::ModelTracker(MilpModel& model)
    : model_(model), deleted_cons_(model.GetNCons()), deleted_vars_(model.GetNVars()) {
    orig_n_vars_ = model.GetNVars();
    orig_var_idx_.resize(orig_n_vars_);
    for (Index iv = 0; iv < orig_n_vars_; ++iv) {
        orig_var_idx_[iv] = iv;
    }

    CalcLocks();
    CalcActivities();
}

void ModelTracker::CompressCons() {
    auto m = model_.GetNCons();
    auto n = model_.GetNVars();

    // Ar & Rhs & Activities
    auto& rhs = model_.GetRhs();

    Index i_write = 0;
    std::vector<Index> new_index_map(m, -1);
    for (Index i_read = 0; i_read < m; ++i_read) {
        if (!deleted_cons_.Get(i_read)) {
            new_index_map[i_read] = i_write;

            if (i_write != i_read) {
                model_.GetRow(i_write) = std::move(model_.GetRow(i_read));
                rhs[i_write] = std::move(rhs[i_read]);
                activities_[i_write] = std::move(activities_[i_read]);
            }
            ++i_write;
        }
    }

    // Ac
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

    model_.Resize(m - deleted_cons_.GetNValues(), n);
    activities_.resize(m - deleted_cons_.GetNValues());

    deleted_cons_.Clear();
}

void ModelTracker::CompressVars() {
    auto m = model_.GetNCons();
    auto n = model_.GetNVars();

    // Objective & Domain & Ac
    auto& coeffs = model_.GetObj().coefficients;
    auto& domain = model_.GetDomain();

    Index i_write = 0;
    std::vector<Index> new_index_map(n, -1);
    for (Index i_read = 0; i_read < n; ++i_read) {
        if (!VarIsDeleted(i_read)) {
            if (i_write != i_read) {
                orig_var_idx_[i_write] = std::move(orig_var_idx_[i_read]);
                coeffs[i_write] = std::move(coeffs[i_read]);
                model_.GetCol(i_write) = std::move(model_.GetCol(i_read));
                domain.Move(i_read, i_write);
                locks_.Move(i_read, i_write);
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
        if (VarIsDeleted(iv1) or VarIsDeleted(iv2)) continue;
        implications_[i].x_ind = new_index_map[iv1];
        implications_[i].y_ind = new_index_map[iv2];
        if (new_n_impl != i) std::swap(implications_[i], implications_[new_n_impl]);
        new_n_impl++;
    }
    implications_.resize(new_n_impl);

    model_.Resize(m, n - deleted_vars_.GetNValues());

    deleted_vars_.Clear();
}

void ModelTracker::CalcLocks() {
    Index n = model_.GetNVars();
    locks_.Resize(n);
    for (Index iv = 0; iv < n; iv++) {
        UpdLocks(iv);
    }
}

void ModelTracker::UpdLocks(Index iv) {
    Index n_up_locks = 0, n_down_locks = 0;
    for (SvIterator el(model_.GetCol(iv)); el; ++el) {
        if (el.value() > 0) {
            if (model_.GetRhs(el.index()).le != -kInf) n_down_locks++;
            if (model_.GetRhs(el.index()).ri != kInf) n_up_locks++;
        }
        if (el.value() < 0) {
            if (model_.GetRhs(el.index()).le != -kInf) n_up_locks++;
            if (model_.GetRhs(el.index()).ri != kInf) n_down_locks++;
        }
    }
    locks_.n_locks[iv][LockType2Index(LockType::kDown)] = n_down_locks;
    locks_.n_locks[iv][LockType2Index(LockType::kUp)] = n_up_locks;
    locks_.sum_locks[iv] = n_down_locks + n_up_locks;
}

void ModelTracker::CalcActivities() {
    auto m = model_.GetNCons();
    activities_.resize(m);
    for (Index ic = 0; ic < m; ic++) {
        UpdActivity(ic);
    }
}

void ModelTracker::UpdActivity(Index ic) {
    activities_[ic] = Activity();
    for (SvIterator el(model_.GetRow(ic)); el; ++el) {
        activities_[ic].AddTerm(el.value(), model_.GetBounds(el.index()));
    }
}

Bounds ModelTracker::DeriveBounds(Index ic, Index iv, Activity act, const Bounds& bnd,
                                  Scalar val) const {
    // Todo: pass non-const bounds and modify them
    // We will modify the activity, so it is passed as a value
    Bounds derived;
    const Bounds& rhs = model_.GetRhs(ic);
    act.RmTerm(val, bnd);
    auto lhs = act.GetRange();
    if (StrongGt(val, 0.0)) {
        derived.le = (rhs.le - lhs.ri) / val;
        derived.ri = (rhs.ri - lhs.le) / val;
    } else if (StrongLt(val, 0.0)) {
        derived.le = (rhs.ri - lhs.le) / val;
        derived.ri = (rhs.le - lhs.ri) / val;
    } else {  // ignore near-zeroes as noise
        return bnd;
    }
    if (model_.GetIntegrality(iv)) {
        derived.le = WeakCeil(derived.le);
        derived.ri = WeakFloor(derived.ri);
    }
    return derived;
}

}  // namespace reshala
