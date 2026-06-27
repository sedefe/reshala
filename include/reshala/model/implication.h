#pragma once

#include "reshala/types.h"

namespace reshala {

struct Implication {
    Index x_ind;
    bool x_value;
    Index y_ind;
    bool left;
    Scalar b;
};

}  // namespace reshala
