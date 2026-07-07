#include "reshala/lp/dual_simplex.h"

namespace reshala {

void DualSimplex::Chuzc() {
    Scalar a_pj, c_j, ratio;
    iv_entering = -1;
    Scalar threshold = kInf;

    for (Index iv = 0; iv < n; iv++) {
        a_pj = a_p[iv] * (s_p * d_n[iv]);
        if (a_pj > kPivotTolerance) {
            c_j = c_n[iv] * d_n[iv];
            ratio = (c_j + kPivotTolerance) / a_pj;
            threshold = std::min(threshold, ratio);
        }
    }

    a_pq = kPivotTolerance;  // Чтобы ниже не сравнивать a_pj с kPivotTolerance
    for (Index iv = 0; iv < n; iv++) {
        a_pj = a_p[iv] * (s_p * d_n[iv]);
        if (a_pj > a_pq) {
            c_j = c_n[iv] * d_n[iv];
            ratio = c_j / a_pj;
            if (ratio <= threshold) {
                a_pq = a_pj;
                c_q = c_j;
                iv_entering = iv;
            }
        }
    }

    if (iv_entering >= 0) {
        theta_p = d_n[iv_entering] * primal_infeasibility / a_pq;
        theta_d = s_p * c_q / a_pq;
    }
}

}  // namespace reshala
