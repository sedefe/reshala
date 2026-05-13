#include "reshala/presolve/transforms.h"

namespace reshala {

void FixVariableTransform::Do(ModelInfo& info) {
    MilpModel& model = info.GetModel();

    model.GetObj().c0 += model.GetObj().coefficients[iv_] * val_;

    for (SvIterator el(model.GetAc().GetCol(iv_)); el; ++el) {
        const Bounds& rhs = model.GetRhs()[el.index()];
        model.GetRhs()[el.index()] = {rhs.le - el.value() * val_, rhs.ri - el.value() * val_};
    }

    info.MaskVar(iv_);
}

void FixVariableTransform::Undo(Solution& sol) { sol.x.insert(sol.x.begin() + iv_, val_); }

}  // namespace reshala
