#include "reshala/milp/branching.h"
#include "reshala/utils.h"

namespace reshala {

Index MostInfeasible::Branch(const Solution& sol) {
    Index candidate = 0;
    Scalar max_fraction = -kInf;

    for (Index iv = 0; iv < sol.x.size(); ++iv) {
        if (!model_.GetIntegrality(iv)) continue;
        Scalar current_fraction = GetFraction(sol.x[iv]);
        if (current_fraction > max_fraction) {
            max_fraction = current_fraction;
            candidate = iv;
        }
    }

    return candidate;
}

}  // namespace reshala
