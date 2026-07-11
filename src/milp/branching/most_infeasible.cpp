#include "reshala/lp/dual_simplex.h"
#include "reshala/milp/branching.h"
#include "reshala/utils.h"

namespace reshala {

Index MostInfeasible::Branch(Node& parent, DualSimplex& ds) {
    Index candidate = -1;
    Scalar max_fraction = -kInf;

    for (Index iv = 0; iv < parent.sol.x.size(); ++iv) {
        if (!model_.GetIntegrality(iv)) continue;
        Scalar current_fraction = GetFraction(parent.sol.x[iv]);
        if (current_fraction > max_fraction) {
            max_fraction = current_fraction;
            candidate = iv;
        }
    }

    Scalar x_floor = Floor(parent.sol.x[candidate]);
    const Bounds orig_bnd = model_.GetBounds(candidate);
    std::array<Bounds, 2> cand_bounds{{{orig_bnd.le, x_floor}, {x_floor + 1, orig_bnd.ri}}};

    Index num_ch = 0;
    for (Index i = 0; i < 2; i++) {
        ds.Restore(parent.ds_state);
        ds.SetBounds(candidate, cand_bounds[i]);
        auto sol = ds.Solve(true);
        children_[i] = Node(parent.level + 1, sol, model_.GetDomain(), ds.Store());

        num_ch += (sol.status == LpStatus::kOptimal);
    }

    return num_ch;
}

}  // namespace reshala
