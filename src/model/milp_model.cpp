#include "reshala/model/milp_model.h"

#include <assert.h>

#include "reshala/utils.h"

namespace reshala {

bool MilpModel::IsIntegerFeasible(const std::vector<Scalar>& x) {
    for (Index iv = 0; iv < GetNVars(); iv++) {
        if (GetIntegrality(iv) and GetFraction(x[iv]) > kEpsZero) {
            return false;
        }
    }
    return true;
}

FeasibilityReport MilpModel::GetFeasReport(const std::vector<Scalar>& x) {
    FeasibilityReport rep;
    rep.max_int_infeas = 0;
    rep.max_bnd_infeas = 0;
    rep.max_con_infeas = 0;

    for (Index iv = 0; iv < GetNVars(); iv++) {
        if (GetIntegrality(iv)) {
            rep.max_int_infeas = std::max(rep.max_int_infeas, GetFraction(x[iv]));
        }
        const Bounds& bnd = GetBounds(iv);
        rep.max_bnd_infeas = std::max(rep.max_bnd_infeas, bnd.le - x[iv]);
        rep.max_bnd_infeas = std::max(rep.max_bnd_infeas, x[iv] - bnd.ri);
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
    for (Index ic = 0; ic < Index(m); ic++) {
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
        return {status, kNan, {}};
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
        if (rhs.le != -kInf) {
            os << lhs << " >= " << rhs.le << std::endl;
        }
        if (rhs.ri != kInf) {
            os << lhs << " <= " << rhs.ri << std::endl;
        }
    }

    os << "Bounds\n";
    for (Index i = 0; i < n; i++) {
        const Bounds& bnd = model.GetBounds(i);
        if (bnd.le != 0.0) {
            os << "x" << i << " >= " << bnd.le << std::endl;
        }
        if (bnd.ri != kInf) {
            os << "x" << i << " <= " << bnd.ri << std::endl;
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
