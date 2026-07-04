#include "reshala/cuts/cut.h"

namespace reshala {

std::ostream& operator<<(std::ostream& os, const Cut& cut) {
    os << cut.lhs << " >= " << cut.rhs << "\n";
    return os;
}

}  // namespace reshala
