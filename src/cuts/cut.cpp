#include "reshala/cuts/cut.h"

namespace reshala {

std::string CutType2Str(CutType type) {
    switch (type) {
        case CutType::kProbing:
            return "Probing";
        case CutType::kCmir:
            return "c-Mir";
        default:
            assert(false && "Unknows cut type");
            return "";
    }
}

std::ostream& operator<<(std::ostream& os, const Cut& cut) {
    os << cut.lhs << " >= " << cut.rhs << "\n";
    return os;
}

}  // namespace reshala
