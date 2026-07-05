#pragma once

#include <vector>

#include "reshala/types.h"

namespace reshala {

struct Scaling {
    std::vector<Index> row;
    std::vector<Index> col;
};

}  // namespace reshala
