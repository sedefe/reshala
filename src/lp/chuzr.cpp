#include "reshala/lp/dual_simplex.h"

namespace reshala {

void DualSimplex::Chuzr() {
    iv_leaving = -1;
    s_p = 0;

    Scalar infeasibility;
    Scalar max_attractiveness = 0.;
    Scalar attractiveness;

    for (Index ic = 0; ic < m; ic++) {
        DenseVector e_p(m, 0.0);
        lina.Btran(SparseVector::UnitVector(m, ic), e_p);
        w_dse[ic] = 0.;
        for (auto x : e_p) w_dse[ic] += x * x;
        w_dse[ic] = std::max(w_dse[ic], kMinDseWeight);
    }

    for (Index ic = 0; ic < m; ic++) {
        const auto& bnd = model_.GetBounds(basis.Basis()[ic]);
        infeasibility = std::max(x_b[ic] - bnd.ri, bnd.le - x_b[ic]);

        if (WeakLe(infeasibility, 0.)) continue;

        attractiveness = infeasibility * infeasibility / w_dse[ic];

        if (attractiveness > max_attractiveness) {
            max_attractiveness = attractiveness;
            iv_leaving = ic;
        }
    }

    if (iv_leaving >= 0) {
        const auto& bnd = model_.GetBounds(basis.Basis()[iv_leaving]);
        if ((x_b[iv_leaving] - bnd.ri) > (bnd.le - x_b[iv_leaving])) {
            s_p = +1;
            primal_infeasibility = x_b[iv_leaving] - bnd.ri;
        } else {
            s_p = -1;
            primal_infeasibility = bnd.le - x_b[iv_leaving];
        }
    }
}

}  // namespace reshala
