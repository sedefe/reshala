#include "reshala/presolve/transforms.h"

namespace reshala {

void FixVariableTransform::Do(MilpModel& model) {
    auto m = model.GetNCons();
    auto n = model.GetNVars();
    // Objective obj_;
    // SparseColMatrix Ac_;
    // SparseRowMatrix Ar_;
    // std::vector<Bounds> rhs_;
    // Domain domain_;

    model.GetObj().coefficients.erase(model.GetObj().coefficients.begin() + iv_);

    for (SparseVector& row : model.GetAr().GetRows()) {
        row.Erase(iv_);
    }

    for (SvIterator el(model.GetAc().GetCol(iv_)); el; ++el) {
        const Bounds& rhs = model.GetRhs()[el.index()];
        model.GetRhs()[el.index()] = {rhs.le - el.value() * val_, rhs.ri - el.value() * val_};
    }

    model.GetAc().GetCols().erase(model.GetAc().GetCols().begin() + iv_);

    model.GetDomain().Erase(iv_);
    model.Resize(m, n - 1);
}

void FixVariableTransform::Undo(Solution& sol) { sol.x.insert(sol.x.begin() + iv_, val_); }

}  // namespace reshala
