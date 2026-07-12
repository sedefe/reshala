#include "reshala/lp/dual_simplex.h"

namespace reshala {

void DualSimplex::SetModel(MilpModel& model) {
    model_orig_ = &model;
    model_ = model;
    m = model_.GetNCons();
    n = model_.GetNVars();
    basis = LpBasis(m, n);
    c_n.resize(n);
    x_b.resize(m);
    e_p.resize(m);
    a_p.resize(n);
    a_q.resize(m);
    d_n.resize(n);

    scaling.ScaleModel(model_);
    lina = Lina(model_.GetAc(), model_.GetAr(), &basis);
}

void DualSimplex::Init() {
    // basic -> non_basis -> c_n -> d_n -> x_b
    basis.Reset();
    lina.Refactor();

    c_n = model_.GetObj().coefficients;

    DenseVector x_n(n, 0);
    for (Index iv = 0; iv < n; iv++) {
        Index i_nb = basis.NonBasis()[iv];
        const Bounds& bnd = model_.GetBounds(i_nb);
        switch (model_.GetType(i_nb)) {
            case BndType::kBoxed:
                if (c_n[iv] >= 0.0) {
                    d_n[iv] = 1;
                    x_n[iv] = bnd.le;
                } else {
                    d_n[iv] = -1;
                    x_n[iv] = bnd.ri;
                }
                break;
            case BndType::kLower:
                d_n[iv] = 1;
                x_n[iv] = bnd.le;
                break;
            case BndType::kUpper:
                d_n[iv] = -1;
                x_n[iv] = bnd.ri;
                break;
            case BndType::kFree:
                d_n[iv] = 1;  // Todo тут ещё a_p[iv] надо смотреть
                x_n[iv] = 0.0;
                break;
            case BndType::kFixed:
                d_n[iv] = 0;
                x_n[iv] = bnd.le;
                break;
            default:
                assert(false);
        }
    }

    MulNRight(x_n, x_b);
    for (Scalar& x : x_b) x = -x;
}

Solution DualSimplex::Solve(bool warm) {
    AddSlacks();
    LpStatus status;

    if (!warm) {
        Init();
    }

    while (true) {
        stats.n_iter += 1;
        // DebugPrint();

        Chuzr();
        if (iv_leaving < 0) {
            status = LpStatus::kOptimal;
            break;
        }
        // std::cout << "Leaving: " << iv_leaving << " (" << basis.Basis()[iv_leaving]
        //           << "), pinf: " << primal_infeasibility << "\n";

        Btran();
        Price();

        Chuzc();
        if (iv_entering < 0) {
            status = LpStatus::kInfeasible;
            break;
        }
        // std::cout << "Entering: " << iv_entering << " (" << basis.NonBasis()[iv_entering]
        //           << "), a_pq: " << a_pq << "\n";

        Ftran();
        Update();
    }

    PruneSlacks();

    DenseVector x;  // Fill & unscale
    if (status == LpStatus::kOptimal) {
        x.resize(n);
        for (Index ic = 0; ic < m; ic++) {
            Index i_b = basis.Basis()[ic];
            if (i_b < n) {
                x[i_b] = std::ldexp(x_b[ic], -scaling.col[i_b]);
            }
        }
        for (Index iv = 0; iv < n; iv++) {
            Index i_nb = basis.NonBasis()[iv];
            if (i_nb < n) {
                x[i_nb] = std::ldexp(GetXnValue(iv), -scaling.col[i_nb]);
            }
        }
    }

    return model_orig_->PrepareSolution(status, x);
}

void DualSimplex::Btran() { lina.Btran(SparseVector(m, iv_leaving, 1.0), e_p); }

void DualSimplex::Price() { MulNLeft(e_p, a_p); }

void DualSimplex::Ftran() { lina.Ftran(model_.GetCol(basis.NonBasis()[iv_entering]), a_q); }

}  // namespace reshala
