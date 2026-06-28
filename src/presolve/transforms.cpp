#include "reshala/presolve/transforms.h"

namespace reshala {

void FixVariableTransform::Undo(Solution& sol) { sol.x[iv_] = val_; }

void SimpleSubTransform::Undo(Solution& sol) { sol.x[iv1_] = a_ * sol.x[iv2_] + b_; }

void LinCombTransform::Undo(Solution& sol) {
    sol.x[iv_] = b_;
    for (SvIterator el(sv_); el; ++el) {
        sol.x[iv_] += el.value() * sol.x[el.index()];
    }
}

}  // namespace reshala
