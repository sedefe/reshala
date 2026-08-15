#include "reshala/lina/basis.h"

namespace reshala {

void LpBasis::AddBasicVars(Index n_vars) {
    Index m = basis.size();
    Index n = non_basis.size();
    basis.resize(m + n_vars);
    index2nb.resize(n + m + n_vars);

    for (Index i = 0; i < n_vars; i++) {
        basis[m + i] = n + m + i;
        index2nb[n + m + i] = -1;
    }
}

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
