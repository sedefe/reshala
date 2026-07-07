#include "reshala/lp/dual_simplex.h"

namespace reshala {

void DualSimplex::Chuzr() {
    iv_leaving = -1;
    s_p = 0;
    Scalar max_infeasibility = kEpsZero;
    Scalar infeasibility;

    for (Index ic = 0; ic < m; ic++) {
        const Bounds& bnd = model_->GetBounds(basis.Basis()[ic]);
        infeasibility = std::max(x_b[ic] - bnd.ri, bnd.le - x_b[ic]);
        if (infeasibility > max_infeasibility) {
            max_infeasibility = infeasibility;
            iv_leaving = ic;
        }
    }

    if (iv_leaving >= 0) {
        auto& bounds = model_->GetBounds(basis.Basis()[iv_leaving]);
        if ((x_b[iv_leaving] - bounds.ri) > (bounds.le - x_b[iv_leaving])) {
            s_p = +1;
            primal_infeasibility = x_b[iv_leaving] - bounds.ri;
        } else {
            s_p = -1;
            primal_infeasibility = bounds.le - x_b[iv_leaving];
        }
    }
}

}  // namespace reshala
