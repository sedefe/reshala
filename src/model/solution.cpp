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
        case LpStatus::kDropped:
            return "Dropped";
        case LpStatus::kUnbounded:
            return "Unbounded";
        case LpStatus::kError:
            return "Error";
        default:
            assert(false && "Unknown status");
            return "";
    }
}

}  // namespace reshala
