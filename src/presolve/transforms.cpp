#include "reshala/presolve/transforms.h"

namespace reshala {

void FixVariableTransform::Undo(Solution& sol) { sol.x[iv_] = val_; }

void SimpleSubTransform::Undo(Solution& sol) { sol.x[iv1_] = a_ * sol.x[iv2_] + b_; }

}  // namespace reshala
