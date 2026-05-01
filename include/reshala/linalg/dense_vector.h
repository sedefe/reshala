#pragma once

#include <ostream>
#include <vector>

#include "reshala/linalg/constants.h"

namespace reshala {
using DenseVector = std::vector<Scalar>;

std::ostream &operator<<(std::ostream &os, const DenseVector &sv);

}  // namespace reshala
