#pragma once

#include "reshala/types.h"

namespace reshala {

struct FeasibilityReport {
    Scalar max_int_infeas;
    Scalar max_bnd_infeas;
    Scalar max_con_infeas;
};

}  // namespace reshala
