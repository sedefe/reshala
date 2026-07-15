#include "reshala/lina/basis.h"

namespace reshala {

std::ostream& operator<<(std::ostream& os, const LpBasis& basis) {
    std::cout << "Basis   : ";
    for (auto ic : basis.Basis()) std::cout << ic << " ";
    std::cout << "\n";
    std::cout << "Nonbasis: ";
    for (auto iv : basis.NonBasis()) std::cout << iv << " ";
    std::cout << "\n";
    return os;
}

}  // namespace reshala
