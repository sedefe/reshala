#include "reshala/presolve/transforms.h"

namespace reshala {

void FixVariableTransform::Undo(Solution& sol) { sol.x[iv_] = val_; }

}  // namespace reshala
