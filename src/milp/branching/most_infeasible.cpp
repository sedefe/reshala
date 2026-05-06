#include "reshala/milp/branching.h"
#include "reshala/utils.h"

namespace reshala {

Index MostInfeasible::Branch(const Solution& sol) {
    Index candidate = 0;
    Scalar max_fraction = GetFraction(sol.x[0]);

    for (Index i = 1; i < sol.x.size(); ++i) {
        Scalar current_fraction = GetFraction(sol.x[i]);
        if (current_fraction > max_fraction) {
            max_fraction = current_fraction;
            candidate = i;
        }
    }

    return candidate;
}

}  // namespace reshala
