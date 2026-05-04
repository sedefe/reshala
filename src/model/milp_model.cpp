#include "reshala/model/milp_model.h"

#include <assert.h>

namespace reshala {

void MilpModel::AddSlacks() {
    assert(!has_slacks_);
    has_slacks_ = true;

    // Adding slacks
    auto m = GetNCons();
    auto n = GetNVars();

    vars_.Resize(n + m);
    for (Index ic = 0; ic < Index(m); ic++) {
        vars_.bounds[n + ic] = {-rhs_[ic].ri, -rhs_[ic].le};
        vars_.types[n + ic] = GetType(vars_.bounds[n + ic]);
        vars_.integrality[n + ic] = false;
    }
}

void MilpModel::PruneSlacks() {
    if (has_slacks_) {
        vars_.Resize(GetNVars());
    }
    has_slacks_ = false;
}

void MilpModel::FinalizeAc() {
    Srm2Scm(Ar_, Ac_);
}

Solution MilpModel::PrepareSolution(const LpStatus status, const std::vector<Scalar> &x) const {
    if (status != LpStatus::kOptimal) {
        return {status, kNan, {}};
    }

    auto res_x = x;
    for (auto &&element : res_x) {
        if (abs(element) < kEpsZero) {
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
        if (bnd.le != -kInf) {
            os << "x[" << i << "] >= " << bnd.le << std::endl;
        }
        if (bnd.ri != kInf) {
            os << "x[" << i << "] <= " << bnd.ri << std::endl;
        }
    }

    os << "Generals\n";
    for (Index i = 0; i < n; i++) {
        if (model.GetVars().integrality[i]) {
            os << i << " ";
        }
    }
    os << "\n";
    os << "End\n";

    return os;
}

}  // namespace reshala
