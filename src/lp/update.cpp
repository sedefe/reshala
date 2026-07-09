#include "reshala/lp/dual_simplex.h"

namespace reshala {

void DualSimplex::RebuildAll() {
    lina.Refactor();

    DenseVector c_b(m);
    DenseVector tmp(m);

    {  // Update c_n
        for (Index ic = 0; ic < m; ic++) {
            Index ib = basis.Basis()[ic];
            c_b[ic] = (ib < n) ? model_->GetObj().coefficients[ib] : 0.0;
        }
        lina.Btran(c_b, tmp);
        MulNLeft(tmp, c_n);
        for (Index iv = 0; iv < n; iv++) {
            Index inb = basis.NonBasis()[iv];
            if (inb < n) {
                c_n[iv] = model_->GetObj().coefficients[inb] - c_n[iv];
            } else {
                c_n[iv] = -c_n[iv];
            }
        }
    }

    {  // Update d_n
        for (Index iv = 0; iv < n; iv++) {
            Index i_nb = basis.NonBasis()[iv];
            const Bounds& bnd = model_->GetBounds(i_nb);
            Scalar d_obj;
            switch (model_->GetType(i_nb)) {
                case BndType::kBoxed:
                    d_obj = (bnd.ri - bnd.le) * c_n[iv];
                    if (StrongGt(c_n[iv], 0.0) or StrongGt(d_obj, 0.0)) {
                        d_n[iv] = 1;
                    } else if (StrongLt(c_n[iv], 0.0) or StrongLt(d_obj, 0.0)) {
                        d_n[iv] = -1;
                    } else {  // Если и кост маленький, и флип ни на что не влияет, оставляем как
                              // есть, чтобы лишний раз не флипать
                    }
                    break;
                case BndType::kLower:
                    d_n[iv] = 1;
                    break;
                case BndType::kUpper:
                    d_n[iv] = -1;
                    break;
                case BndType::kFree:
                    d_n[iv] = 1;  // Todo тут ещё a_p[iv] надо смотреть
                    break;
                case BndType::kFixed:
                    d_n[iv] = 0;
                    break;
                default:
                    assert(false);
            }
        }
        BndType type = model_->GetType(basis.NonBasis()[iv_entering]);
        d_n[iv_entering] = (type == BndType::kFixed) ? 0 : (s_p > 0 ? -1 : 1);
    }

    {  // Updade x_b
        DenseVector x_n(n, 0.0);
        for (Index iv = 0; iv < n; iv++) x_n[iv] = GetXnValue(iv);
        MulNRight(x_n, tmp);
        lina.Ftran(tmp, x_b);
        for (Scalar& x : x_b) x = -x;
    }
}

void DualSimplex::Update() {
    auto x_q_old = GetXnValue(iv_entering);
    basis.Swap(iv_leaving, iv_entering);

    if (lina.GetAge() >= kMaxLinaAge) {
        RebuildAll();
        return;
    }
    lina.Update(iv_leaving, iv_entering);  // Update lina

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
        d_n[iv_entering] = (type == BndType::kFixed) ? 0 : (s_p > 0 ? -1 : 1);
    }

    {  // Update x_b
        for (Index ic = 0; ic < m; ic++) {
            x_b[ic] -= theta_p * a_q[ic];
        }
        x_b[iv_leaving] = theta_p + x_q_old;
    }
}

}  // namespace reshala
