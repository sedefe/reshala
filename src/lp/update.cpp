#include "reshala/lp/dual_simplex.h"

namespace reshala {

void DualSimplex::Update() {
    auto x_q_old = GetXnValue(iv_entering);

    {  // Update lina
        basis.Swap(iv_leaving, iv_entering);

        if (lina.GetAge() < kMaxLinaAge) {
            lina.Update(iv_leaving, iv_entering);
        } else {
            lina.Refactor();
        }
    }

    {  // Update c_n
        for (Index iv = 0; iv < n; iv++) {
            Scalar old = c_n[iv];
            c_n[iv] -= theta_d * a_p[iv];
            if (old * c_n[iv] < -0.1 and iv != iv_entering and
                model_->GetType(basis.NonBasis()[iv]) != BndType::kFixed) {
                // Might make us dual infeasible
                std::cerr << "Abnormal c_n update: " << old << " -> " << c_n[iv] << "\n";
            }
        }
        c_n[iv_entering] = -theta_d;
    }

    {  // Update d_n
        const Bounds& bnd = model_->GetBounds(basis.NonBasis()[iv_entering]);
        BndType type = model_->GetType(basis.NonBasis()[iv_entering]);
        if (type == BndType::kFixed) {
            d_n[iv_entering] = 0;  // A fixed can't enter the basis
        } else {
            if (s_p > 0) {
                d_n[iv_entering] = -1;
            } else {
                d_n[iv_entering] = 1;
            }
        }
    }

    {  // Update x_b
        for (Index ic = 0; ic < m; ic++) {
            x_b[ic] -= theta_p * a_q[ic];
        }
        x_b[iv_leaving] = theta_p + x_q_old;
    }
}

}  // namespace reshala
