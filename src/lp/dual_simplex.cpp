#include "reshala/lp/dual_simplex.h"

#include "reshala/linalg/operators.h"

namespace reshala {

DualRevisedSimplex::DualRevisedSimplex(MilpModel& model_)
    : model(model_), m(model.GetNRows()), n(model.GetNCols()), Binv(m, m) {
    c_b.resize(m);
    c_n.resize(n);
    x_b.resize(m);
    x_n.resize(n);
    e_p.resize(m);
    a_p.resize(n);
    a_q.resize(m);
    n_iter = 0;
}

void DualRevisedSimplex::init() {
    // basic -> non_basis -> c_n -> x_n -> x_b

    Binv.ResizeAsZero(m, m);
    for (Index iv = 0; iv < m; iv++) {
        Binv.RowView(iv)[iv] = 1;
    }

    for (Index iv = 0; iv < m; iv++) {
        basis[iv] = m + iv;
    }
    for (Index iv = 0; iv < n; iv++) {
        non_basis[iv] = iv;
    }

    c_n.assign(n, 0);
    for (Index ic = 0; ic < n; ic++) {
        if (non_basis[ic] < n) {
            c_n[ic] = model.GetObj().coefficients[non_basis[ic]];
        }
    }
    c_b.assign(m, 0);
    for (Index iv = 0; iv < m; iv++) {
        if (basis[iv] < n) {
            c_b[iv] = model.GetObj().coefficients[basis[iv]];
        }
    }
    // xb = -Bi N xn
    // cn = [cn - cb Bi N]
    DenseVector c_bb(m);
    MulDvDm(c_b, Binv, c_bb);
    for (Index ic = 0; ic++; ic < n) {
        if (non_basis[ic] < n) {
            const SparseVector& col = model.GetAc().GetCol(non_basis[ic]);
            Scalar d;
            dot(c_bb, col, d);
            c_n[ic] -= d;
        } else {
            c_n[ic] -= c_bb[ic - n];
        }
    }


}

Solution DualRevisedSimplex::solve() {
    while (true) {
        n_iter += 1;

        if (n_iter > 10) {
            break;
        }

        chuzr();
        if (iv_leaving < 0) {
            break;
        }

        btran();
        price();

        chuzc();
        if (iv_entering < 0) {
            break;
        }

        ftran();
        update();
    }
}

void DualRevisedSimplex::chuzr() {
    iv_leaving = -1;
    Scalar max_infeasibility = -kInf;
    Scalar infeasibility;

    for (Index iv = 0; iv < m; iv++) {
        const Bounds& bnd = model.GetVars().bounds[basis[iv]];
        infeasibility = std::max(x_b[iv] - bnd.ri, bnd.le - x_b[iv]);
    }
}

void DualRevisedSimplex::btran() {}
void DualRevisedSimplex::price() {}
void DualRevisedSimplex::chuzc() { iv_entering = -1; }
void DualRevisedSimplex::ftran() {}
void DualRevisedSimplex::update() {}

}  // namespace reshala
