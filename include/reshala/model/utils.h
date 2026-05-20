#pragma once

#include "reshala/types.h"

namespace reshala {

struct FeasibilityReport {
    Scalar abs_int_infeas = 0;
    Scalar abs_bnd_infeas = 0;
    Scalar abs_con_infeas = 0;

    Scalar rel_bnd_infeas = 0;
    Scalar rel_con_infeas = 0;
};

std::ostream& operator<<(std::ostream& os, const FeasibilityReport& rep);

}  // namespace reshala
