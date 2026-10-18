#include "reshala/presolve/model_tracker.h"

namespace reshala {

void ModelTracker::FixVar(Index iv, Scalar val) {
    model_.GetObj().c0 += model_.GetObj().coefficients[iv] * val;

    const Bounds& old_bnd = model_.GetBounds(iv);
    for (SvIterator el(model_.GetCol(iv)); el; ++el) {
        Index ic = el.index();
        activities_[ic].RmTerm(el.value(), old_bnd);

        const Bounds& rhs = model_.GetRhs(ic);
        model_.GetRhs(ic) = {rhs.le - el.value() * val, rhs.ri - el.value() * val};
    }

    transforms_.push_back(std::make_unique<FixVariableTransform>(orig_var_idx_[iv], val));
    MaskVar(iv);
}

void ModelTracker::ConstShiftVar(Index iv, Scalar val) {
    model_.GetObj().c0 -= model_.GetObj().coefficients[iv] * val;

    const Bounds& bnd = model_.GetBounds(iv);
    UpdVarBounds(iv, {bnd.le + val, bnd.ri + val});

    for (SvIterator el(model_.GetCol(iv)); el; ++el) {
        const Bounds& rhs = model_.GetRhs(el.index());
        model_.GetRhs(el.index()) = {rhs.le + el.value() * val, rhs.ri + el.value() * val};
    }

    transforms_.push_back(std::make_unique<ConstShiftTransform>(orig_var_idx_[iv], val));
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
        if (ConIsDeleted(ic)) continue;

        Scalar val_iv1 = el.value();
        Scalar val_iv2 = a * val_iv1;
        SparseVector& row = model_.GetRow(ic);
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
        const Bounds& rhs = model_.GetRhs(ic);
        model_.GetRhs(ic) = {rhs.le - val_iv1 * b, rhs.ri - val_iv1 * b};
    }
    // Ac
    model_.GetCol(iv2) = axpy(a, model_.GetCol(iv1), model_.GetCol(iv2));

    if (StrongGt(new_bnd2.le, bnd2.le) or StrongLt(new_bnd2.ri, bnd2.ri)) {
        UpdVarBounds(iv2, std::move(new_bnd2));
    }

    UpdLocks(iv2);  // Todo-locks: can update in the loop above

    transforms_.push_back(
        std::make_unique<SimpleSubTransform>(orig_var_idx_[iv1], a, orig_var_idx_[iv2], b));
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
        if (VarIsDeleted(el.index())) continue;
        model_.GetObj().coefficients[el.index()] -= c * el.value() / a;
        sv.Push(orig_var_idx_[el.index()], -el.value() / a);
    }

    UpdVarBounds(iv, {0., 0.});

    transforms_.push_back(std::make_unique<LinCombTransform>(orig_var_idx_[iv], sv, b / a));
    MaskVar(iv);
}

void ModelTracker::UpdCon(Index ic, const SparseVector& row) { model_.GetRow(ic) = row; }

void ModelTracker::RebuildAc() { model_.FinalizeAc(); }

void ModelTracker::UpdRhs(Index ic, Bounds rhs) {
    if (WeakEq(rhs.le, rhs.ri)) rhs.le = rhs.ri = (rhs.le + rhs.ri) / 2;

    model_.GetRhs(ic) = rhs;

    for (SvIterator el(model_.GetRow(ic)); el; ++el) {
        UpdLocks(el.index());  // Todo-locks: check old&new rhs
    }

    stat.n_ch_rhs++;
}

void ModelTracker::UpdVarBounds(Index iv, Bounds bnd) {
    const Bounds& old_bnd = model_.GetBounds(iv);
    if (WeakEq(bnd.le, bnd.ri)) bnd.le = bnd.ri = (bnd.le + bnd.ri) / 2;

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
    if (IsZero(val)) {
        model_.GetRow(ic).EraseIndex(iv);
        model_.GetCol(iv).EraseIndex(ic);
        UpdLocks(iv);  // Todo-locks: look at val
    }
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

void ModelTracker::ScaleCon(Index ic, Scalar scale) {
    SparseVector& row = model_.GetRow(ic);
    row *= scale;
    for (SvIterator el(row); el; ++el) {
        auto iv = el.index();
        if (VarIsDeleted(iv)) continue;
        model_.GetCol(iv).AtRef(ic) *= scale;
    }
    activities_[ic].Scale(scale);

    const auto& rhs = model_.GetRhs(ic);
    UpdRhs(ic, {rhs.le * scale, rhs.ri * scale});

    stat.n_ch_coeff += row.Size();
    stat.n_ch_rhs++;
}

void ModelTracker::ScaleVar(Index iv, Scalar scale) {
    assert(StrongGt(scale, 0.0));
    SparseVector& col = model_.GetCol(iv);
    col *= scale;
    for (SvIterator el(col); el; ++el) {
        auto ic = el.index();
        if (ConIsDeleted(ic)) continue;
        model_.GetRow(ic).AtRef(iv) *= scale;
    }

    const auto& bnd = model_.GetBounds(iv);
    model_.SetBounds(iv, {bnd.le / scale, bnd.ri / scale});

    model_.GetObj().coefficients[iv] *= scale;

    stat.n_ch_coeff += col.Size();
    stat.n_ch_bnd++;

    transforms_.push_back(std::make_unique<ScaleTransform>(orig_var_idx_[iv], scale));
}

void ModelTracker::SetVarInt(Index iv) {
    const Bounds& bnd = model_.GetBounds(iv);
    model_.SetIntegrality(iv, true);
    model_.SetBounds(iv, {WeakCeil(bnd.le), WeakFloor(bnd.ri)});
}

void ModelTracker::ImportBounder(Bounder& bounder) {
    std::swap(activities_, bounder.activities);
    std::swap(model_.GetDomain(), bounder.domain);
    for (auto ic : bounder.redundant_con_mask.GetValues()) {  // Mini 3.1
        MaskCon(ic);
    }
    for (auto iv : bounder.all_changed_vars.GetValues()) {  // Mini 4.1
        if (model_.GetType(iv) == BndType::kFixed) {
            FixVar(iv, (model_.GetBounds(iv).le + model_.GetBounds(iv).ri) / 2);
        }
    }
    CalcLocks();  // Todo-locks: only for changed cons' vars
}

}  // namespace reshala
