#include "reshala/model/solution.h"

#include <assert.h>

namespace reshala {

std::string LpStatus2Str(LpStatus status) {
    switch (status) {
        case LpStatus::kUnknown:
            return "Unknown";
        case LpStatus::kOptimal:
            return "Optimal";
        case LpStatus::kInfeasible:
            return "Infeasible";
        case LpStatus::kUnbounded:
            return "Unbounded";
        default:
            assert(false && "Unknows status");
    }
}

}  // namespace reshala
