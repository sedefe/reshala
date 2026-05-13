#include "reshala/presolve/transforms.h"

namespace reshala {

void FixVariableTransform::Undo(Solution& sol) { sol.x.insert(sol.x.begin() + iv_, val_); }

}  // namespace reshala
