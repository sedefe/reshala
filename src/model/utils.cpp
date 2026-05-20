#include "reshala/model/utils.h"

namespace reshala {

std::ostream& operator<<(std::ostream& os, const FeasibilityReport& rep) {
    os << "Absolute infeasibilities:\n"
       << "  Integrality: " << rep.abs_int_infeas << "\n"
       << "  Bounds:      " << rep.abs_bnd_infeas << "\n"
       << "  Constraints: " << rep.abs_con_infeas << "\n"
       << "Relative infeasibilities:\n"
       << "  Bounds:      " << rep.rel_bnd_infeas << "\n"
       << "  Constraints: " << rep.rel_con_infeas;
    return os;
}

}  // namespace reshala
