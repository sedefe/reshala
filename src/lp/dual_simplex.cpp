#include "reshala/lp/dual_simplex.h"

namespace reshala {

std::ostream& operator<<(std::ostream& os, const DsStats& stats) {
    os << "Lp iters: " << stats.n_iter << "\n";
    return os;
}

void DualSimplex::SetModel(MilpModel& model) {
    model_ = &model;
    m = model_->GetNCons();
    n = model_->GetNVars();
    basis = LpBasis(m, n);
    lina = Lina(&model.GetAc(), &model.GetAr(), &basis);
    c_n.resize(n);
    x_b.resize(m);
    e_p.resize(m);
    a_p.resize(n);
    a_q.resize(m);
    d_n.resize(n);
}

DsState DualSimplex::Store() const { return {c_n, x_b, d_n, basis, lina}; }

void DualSimplex::Restore(const DsState& state) {
    c_n = state.c_n;
    x_b = state.x_b;
    d_n = state.d_n;
    basis = state.basis;
    lina = state.lina;
}

void DualSimplex::Init() {
    // basic -> non_basis -> c_n -> d_n -> x_b
    basis.Reset();
    lina.Reset();

    c_n = model_->GetObj().coefficients;

    DenseVector x_n(n, 0);
    for (Index iv = 0; iv < n; iv++) {
        Index i_nb = basis.NonBasis()[iv];
        const Bounds& bnd = model_->GetBounds(i_nb);
        switch (model_->GetType(i_nb)) {
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

    MulScmDv(model_->GetAc(), x_n, x_b);
    for (Scalar& x : x_b) x = -x;
}

Solution DualSimplex::Solve(bool warm) {
    model_->AddSlacks();
    ForceBounds();
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
        // std::cout << "Entering: " << iv_entering << " (" << basis.NonBasis()[iv_entering] <<
        // ")\n";

        Ftran();
        Update();
    }

    model_->PruneSlacks();
    UnforceBounds();

    DenseVector x;
    if (status == LpStatus::kOptimal) {
        x.resize(n);
        for (Index ic = 0; ic < m; ic++) {
            Index i_b = basis.Basis()[ic];
            if (i_b < n) {
                x[i_b] = x_b[ic];
            }
        }
        for (Index iv = 0; iv < n; iv++) {
            Index i_nb = basis.NonBasis()[iv];
            if (i_nb < n) {
                x[i_nb] = GetXnValue(iv);
            }
        }
    }

    return model_->PrepareSolution(status, x);
}

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

void DualSimplex::Btran() { lina.Btran(iv_leaving, e_p); }

void DualSimplex::Price() {
    a_p.assign(n, 0.0);
    for (Index ic = 0; ic < m; ic++) {
        if (IsZero(e_p[ic])) continue;
        for (SvIterator el(model_->GetRow(ic)); el; ++el) {
            if (basis.Index2Nb()[el.index()] >= 0) {
                a_p[basis.Index2Nb()[el.index()]] += e_p[ic] * el.value();
            }
        }
        if (basis.Index2Nb()[n + ic] >= 0) {
            a_p[basis.Index2Nb()[n + ic]] += e_p[ic];
        }
    }
}

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

void DualSimplex::Ftran() { lina.Ftran(basis.NonBasis()[iv_entering], a_q); }

void DualSimplex::Update() {
    auto x_q_old = GetXnValue(iv_entering);

    {  // Update lina
        basis.Swap(iv_leaving, iv_entering);
        lina.Update(iv_leaving, iv_entering);
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

inline Scalar DualSimplex::GetXnValue(Index iv) {
    // Todo обрабатывать свободные переменные, глядя на a_p
    const Bounds& bnd = model_->GetBounds(basis.NonBasis()[iv]);
    return (d_n[iv] >= 0) ? bnd.le : bnd.ri;
}

void DualSimplex::ForceBounds() {
    initial_domain = model_->GetDomain();
    for (Index iv = 0; iv < m + n; iv++) {
        model_->SetBounds(iv, BoundsIntersection(model_->GetBounds(iv), {-kMaxAbs, kMaxAbs}));
    }
}

void DualSimplex::UnforceBounds() { model_->SetDomain(initial_domain); }

void DualSimplex::DebugPrint() {
    DenseVector x(n);
    for (Index iv = 0; iv < m; iv++) {
        if (basis.Basis()[iv] < n) x[basis.Basis()[iv]] = x_b[iv];
    }
    for (Index iv = 0; iv < n; iv++) {
        if (basis.NonBasis()[iv] < n) x[basis.NonBasis()[iv]] = GetXnValue(iv);
    }
    auto res = model_->PrepareSolution(LpStatus::kOptimal, x);

    std::cout << "===== " << stats.n_iter << " y=" << res.y << " =====\n";
    std::cout << "Basis   : ";
    for (auto ic : basis.Basis()) std::cout << ic << " ";
    std::cout << "\n";
    std::cout << "Nonbasis: ";
    for (auto iv : basis.NonBasis()) std::cout << iv << " ";
    std::cout << "\n";
    std::cout << "c_n: ";
    for (auto x : c_n) std::cout << x << " ";
    std::cout << "\n";
    std::cout << "d_n: ";
    for (auto x : d_n) std::cout << Index(x) << " ";
    std::cout << "\n";
    std::cout << "x_b: ";
    for (auto x : x_b) std::cout << x << " ";
    std::cout << "\n";
}

}  // namespace reshala
