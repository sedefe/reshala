#include "reshala/model/milp_model.h"

#include <assert.h>

#include "reshala/utils.h"

namespace reshala {

Index MilpModel::GetNnz() const {
    Index res = 0;
    for (const auto& row : Ar_.GetRows()) res += row.Size();
    return res;
}

bool MilpModel::IsIntegerFeasible(const std::vector<Scalar>& x) const {
    assert(x.size() == GetNVars());
    for (Index iv = 0; iv < GetNVars(); iv++) {
        if (GetIntegrality(iv) and GetFraction(x[iv]) > kEpsZero) {
            return false;
        }
    }
    return true;
}

FeasibilityReport MilpModel::GetFeasReport(const std::vector<Scalar>& x) const {
    FeasibilityReport rep;
    Scalar abs_infeas, rel_infeas;

    auto CheckBounds = [](const Bounds& bnd, Scalar value, Scalar& max_abs_infeas,
                          Scalar& max_rel_infeas) {
        Scalar le_inf = std::max(bnd.le - value, 0.0);
        Scalar ri_inf = std::max(value - bnd.ri, 0.0);

        Scalar abs_infeas = std::max(le_inf, ri_inf);
        Scalar rel_infeas = std::max(le_inf / std::max(std::abs(bnd.le), 1.0),
                                     ri_inf / std::max(std::abs(bnd.ri), 1.0));

        max_abs_infeas = std::max(max_abs_infeas, abs_infeas);
        max_rel_infeas = std::max(max_rel_infeas, rel_infeas);
    };

    for (Index iv = 0; iv < GetNVars(); iv++) {
        if (GetIntegrality(iv)) {
            rep.abs_int_infeas = std::max(rep.abs_int_infeas, GetFraction(x[iv]));
        }
        CheckBounds(GetBounds(iv), x[iv], rep.abs_bnd_infeas, rep.rel_bnd_infeas);
    }

    auto m = GetNCons();
    DenseVector y(m);
    MulScmDv(Ac_, x, y);
    for (Index ic = 0; ic < m; ic++) {
        CheckBounds(rhs_[ic], y[ic], rep.abs_con_infeas, rep.rel_con_infeas);
    }

    return rep;
}

void MilpModel::AddSlacks() {
    assert(!has_slacks_);
    has_slacks_ = true;

    // Adding slacks
    auto m = GetNCons();
    auto n = GetNVars();

    domain_.Resize(n + m);
    for (Index ic = 0; ic < m; ic++) {
        domain_.SetBounds(n + ic, {-rhs_[ic].ri, -rhs_[ic].le});
        domain_.SetIntegrality(n + ic, false);
    }
}

void MilpModel::PruneSlacks() {
    if (has_slacks_) {
        domain_.Resize(GetNVars());
    }
    has_slacks_ = false;
}

void MilpModel::FinalizeAc() { Srm2Scm(Ar_, Ac_); }

Solution MilpModel::PrepareSolution(const LpStatus status, const std::vector<Scalar>& x) const {
    if (status != LpStatus::kOptimal) {
        return {status, kInf, {}};
    }

    auto res_x = x;
    for (auto&& element : res_x) {
        if (IsZero(element)) {
            element = 0.0;
        }
    }

    return {status, obj_.evaluate(res_x), res_x};
}

std::ostream& operator<<(std::ostream& os, const MilpModel& model) {
    os << model.obj_ << std::endl;

    os << "Subject to\n";
    Index m = model.GetNCons();
    Index n = model.GetNVars();
    for (Index i = 0; i < m; i++) {
        const SparseVector& lhs = model.Ar_.GetRows()[i];
        const Bounds& rhs = model.rhs_[i];

        if (rhs.le == rhs.ri) {
            os << lhs << " = " << rhs.le << std::endl;
        } else {
            if (rhs.le != -kInf) os << lhs << " >= " << rhs.le << std::endl;
            if (rhs.ri != kInf) os << lhs << " <= " << rhs.ri << std::endl;
        }
    }

    os << "Bounds\n";
    for (Index i = 0; i < n; i++) {
        const Bounds& bnd = model.GetBounds(i);
        if (model.GetType(i) == BndType::kFixed) {
            os << "x" << i << " = " << (bnd.le + bnd.ri) / 2 << std::endl;
        } else {
            if (bnd.le != 0.0) os << "x" << i << " >= " << bnd.le << std::endl;
            if (bnd.ri != kInf) os << "x" << i << " <= " << bnd.ri << std::endl;
        }
    }

    os << "Generals\n";
    for (Index i = 0; i < n; i++) {
        if (model.GetIntegrality(i)) {
            os << "x" << i << " ";
        }
    }
    os << "\n";
    os << "End\n";

    return os;
}

}  // namespace reshala
