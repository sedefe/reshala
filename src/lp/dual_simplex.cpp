#include "reshala/lp/dual_simplex.h"

#include "reshala/linalg/operators.h"

namespace reshala {

DualSimplex::DualSimplex(MilpModel& model)
    : model_(model), m(model_.GetNCons()), n(model_.GetNVars()), Binv(m, m) {
    basis.resize(m);
    non_basis.resize(m);
    index2nb.resize(m + n);
    c_n.resize(n);
    x_b.resize(m);
    x_n.resize(n);
    e_p.resize(m);
    a_p.resize(n);
    a_q.resize(m);
    d_p.resize(n);
    n_iter = 0;
}

void DualSimplex::Init() {
    // basic -> non_basis -> c_n -> x_n -> x_b
    basis.resize(m);
    non_basis.resize(n);
    index2nb.resize(n + m);
    for (Index ic = 0; ic < m; ic++) {
        basis[ic] = n + ic;
        index2nb[n + ic] = -1;
    }
    for (Index iv = 0; iv < n; iv++) {
        non_basis[iv] = iv;
        index2nb[iv] = iv;
    }

    Binv.ResizeAsZero(m, m);
    for (Index iv = 0; iv < m; iv++) {
        Binv.RowView(iv)[iv] = 1;
    }

    c_n = model_.GetObj().coefficients;

    x_n = std::vector<Scalar>(n, 0);
    for (Index iv = 0; iv < n; iv++) {
        SetXnValue(iv);
    }

    MulScmDv(model_.GetAc(), x_n, x_b);
    for (Scalar& x : x_b) x = -x;
}

Solution DualSimplex::Solve() {
    model_.AddSlacks();
    ForceBounds();
    LpStatus status;

    Init();
    while (true) {
        n_iter += 1;
        // DebugPrint();

        Chuzr();
        if (iv_leaving < 0) {
            status = LpStatus::kOptimal;
            break;
        }
        // std::cout << "Leaving: " << iv_leaving << " (" << basis[iv_leaving]
        //           << "), pinf: " << primal_infeasibility << "\n";

        Btran();
        Price();

        Chuzc();
        if (iv_entering < 0) {
            status = LpStatus::kInfeasible;
            break;
        }
        // std::cout << "Entering: " << iv_entering << " (" << non_basis[iv_entering] << "\n";

        Ftran();
        Update();
    }

    model_.PruneSlacks();
    UnforceBounds();

    DenseVector x;
    if (status == LpStatus::kOptimal) {
        x.resize(n);
        for (Index ic = 0; ic < m; ic++) {
            if (basis[ic] < n) {
                x[basis[ic]] = x_b[ic];
            }
        }
        for (Index iv = 0; iv < n; iv++) {
            if (non_basis[iv] < n) {
                x[non_basis[iv]] = x_n[iv];
            }
        }
    }

    return model_.PrepareSolution(status, x);
}

void DualSimplex::Chuzr() {
    iv_leaving = -1;
    s_p = 0;
    Scalar max_infeasibility = kEpsZero;
    Scalar infeasibility;

    for (Index ic = 0; ic < m; ic++) {
        const Bounds& bnd = model_.GetBounds(basis[ic]);
        infeasibility = std::max(x_b[ic] - bnd.ri, bnd.le - x_b[ic]);
        if (infeasibility > max_infeasibility) {
            max_infeasibility = infeasibility;
            iv_leaving = ic;
        }
    }

    if (iv_leaving >= 0) {
        auto& bounds = model_.GetBounds(basis[iv_leaving]);
        if ((x_b[iv_leaving] - bounds.ri) > (bounds.le - x_b[iv_leaving])) {
            s_p = +1;
            primal_infeasibility = x_b[iv_leaving] - bounds.ri;
        } else {
            s_p = -1;
            primal_infeasibility = bounds.le - x_b[iv_leaving];
        }
    }
}

void DualSimplex::Btran() { e_p.assign(Binv.RowView(iv_leaving), Binv.RowView(iv_leaving) + m); }

void DualSimplex::Price() {
    // Todo: row-wise pricing
    for (Index iv = 0; iv < n; iv++) {
        if (non_basis[iv] < n) {
            dot(e_p, model_.GetAc().GetCol(non_basis[iv]), a_p[iv]);
        } else {
            a_p[iv] = e_p[non_basis[iv] - n];
        }
    }
}

void DualSimplex::Chuzc() {
    Scalar a_pj, c_j;
    Scalar min_ratio = kInf;
    iv_entering = -1;

    for (Index iv = 0; iv < n; iv++) {
        if (model_.GetType(non_basis[iv]) == VarType::kFixed) {
            continue;
        }
        a_pj = a_p[iv] * (s_p * d_p[iv]);
        c_j = c_n[iv] * d_p[iv];

        if (a_pj > kEpsZero) {
            auto ratio = c_j / a_pj;
            if (ratio < min_ratio) {
                min_ratio = ratio;
                iv_entering = iv;
                a_pq = a_pj;
                c_j_entering = c_j;
            }
        }
    }

    if (iv_entering >= 0) {
        theta_p = d_p[iv_entering] * primal_infeasibility / a_pq;
        theta_d = s_p * c_j_entering / a_pq;
    }
}

void DualSimplex::Ftran() {
    if (non_basis[iv_entering] < n) {
        MulDmSv(Binv, model_.GetAc().GetCol(non_basis[iv_entering]), a_q);
    } else {
        MulDmSv(Binv, SparseVector(m, non_basis[iv_entering] - n, 1.0), a_q);
    }
}

void DualSimplex::Update() {
    DenseVector diff(m, 0.0);
    DenseVector B_p(m, 0);

    Index ib = basis[iv_leaving];
    Index inb = non_basis[iv_entering];

    const SparseVector& leaving_col =
        ib < n ? model_.GetAc().GetCol(ib) : SparseVector(m, ib - n, 1.0);
    const SparseVector& entering_col =
        inb < n ? model_.GetAc().GetCol(inb) : SparseVector(m, inb - n, 1.0);

    {  // Update Binv
        const SparseVector delta = sub(entering_col, leaving_col);
        const DenseVector row(Binv.RowView(iv_leaving), Binv.RowView(iv_leaving) + m);
        Scalar multiplier;
        dot(row, delta, multiplier);
        multiplier = 1 / (1 + multiplier);

        DenseVector d(m, 0);
        MulDmSv(Binv, delta, d);

        for (Index i = 0; i < m; i++) {
            for (Index j = 0; j < m; j++) {
                Binv.RowView(i)[j] -= d[i] * row[j] * multiplier;
            }
        }
    }
    std::swap(basis[iv_leaving], non_basis[iv_entering]);

    {  // Update c_n
        for (Index iv = 0; iv < n; iv++) {
            c_n[iv] -= theta_d * a_p[iv];
        }
        c_n[iv_entering] = -theta_d;
    }

    auto x_q_old = x_n[iv_entering];
    {  // Update x_n
        const Bounds& bnd = model_.GetBounds(non_basis[iv_entering]);
        if (s_p > 0) {
            d_p[iv_entering] = -1;
            x_n[iv_entering] = bnd.ri;
        } else {
            d_p[iv_entering] = 1;
            x_n[iv_entering] = bnd.le;
        }
    }

    {  // Update x_b
        for (Index ic = 0; ic < m; ic++) {
            x_b[ic] -= theta_p * a_q[ic];
        }
        x_b[iv_leaving] = theta_p + x_q_old;
    }
}

void DualSimplex::ForceBounds() {
    initial_domain = model_.GetDomain();
    for (Index iv = 0; iv < m + n; iv++) {
        model_.SetBounds(iv, BoundsIntersection(model_.GetBounds(iv), {-kMaxAbs, kMaxAbs}));
    }
}

void DualSimplex::UnforceBounds() { model_.SetDomain(initial_domain); }

void DualSimplex::SetXnValue(Index iv) {
    const Bounds& bnd = model_.GetBounds(non_basis[iv]);
    switch (model_.GetType(non_basis[iv])) {
        case VarType::kBoxed:
            if (c_n[iv] >= 0.0) {
                d_p[iv] = 1;
                x_n[iv] = bnd.le;
            } else {
                d_p[iv] = -1;
                x_n[iv] = bnd.ri;
            }
            break;
        case VarType::kLower:
            d_p[iv] = 1;
            x_n[iv] = bnd.le;
            break;
        case VarType::kUpper:
            d_p[iv] = -1;
            x_n[iv] = bnd.ri;
            break;
        case VarType::kFree:
            d_p[iv] = 1;  // Todo тут ещё a_p[iv] надо смотреть
            x_n[iv] = 0.0;
            break;
        case VarType::kFixed:
            d_p[iv] = 0;
            x_n[iv] = bnd.le;
            break;
        default:
            assert(false);
    }
}

void DualSimplex::DebugPrint() {
    DenseVector x(n);
    for (Index iv = 0; iv < m; iv++) {
        if (basis[iv] < n) x[basis[iv]] = x_b[iv];
    }
    for (Index iv = 0; iv < n; iv++) {
        if (non_basis[iv] < n) x[non_basis[iv]] = x_n[iv];
    }
    auto res = model_.PrepareSolution(LpStatus::kOptimal, x);

    std::cout << "===== " << n_iter << " y=" << res.y << " =====\n";
    std::cout << "Basis   : ";
    for (auto ic : basis) std::cout << ic << " ";
    std::cout << "\n";
    std::cout << "Nonbasis: ";
    for (auto iv : non_basis) std::cout << iv << " ";
    std::cout << "\n";
    std::cout << "c_n: ";
    for (auto x : c_n) std::cout << x << " ";
    std::cout << "\n";
    std::cout << "x_n: ";
    for (auto x : x_n) std::cout << x << " ";
    std::cout << "\n";
    std::cout << "x_b: ";
    for (auto x : x_b) std::cout << x << " ";
    std::cout << "\n";
}

}  // namespace reshala
