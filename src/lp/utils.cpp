#include "reshala/lp/dual_simplex.h"

namespace reshala {

std::ostream& operator<<(std::ostream& os, const DsStats& stats) {
    os << "Lp iters: " << stats.n_iter << "\n";
    return os;
}

DsState DualSimplex::Store() const { return {c_n, x_b, d_n, basis, lina}; }

void DualSimplex::Restore(const DsState& state) {
    c_n = state.c_n;
    x_b = state.x_b;
    d_n = state.d_n;
    basis = state.basis;
    lina = state.lina;
}

Scalar DualSimplex::GetXnValue(Index iv) {
    // Todo обрабатывать свободные переменные, глядя на a_p
    const Bounds& bnd = model_->GetBounds(basis.NonBasis()[iv]);
    return (d_n[iv] >= 0) ? bnd.le : bnd.ri;
}

void DualSimplex::MulNLeft(const DenseVector& x, DenseVector& res) const {
    res.assign(n, 0.0);
    for (Index ic = 0; ic < m; ic++) {
        if (IsZero(x[ic])) continue;
        for (SvIterator el(model_->GetRow(ic)); el; ++el) {
            if (basis.Index2Nb()[el.index()] >= 0) {
                res[basis.Index2Nb()[el.index()]] += x[ic] * el.value();
            }
        }
        if (basis.Index2Nb()[n + ic] >= 0) {
            res[basis.Index2Nb()[n + ic]] += x[ic];
        }
    }
}

void DualSimplex::MulNRight(const DenseVector& x, DenseVector& res) const {
    res.assign(m, 0.0);
    for (Index iv = 0; iv < n; iv++) {
        if (IsZero(x[iv])) continue;
        Index inb = basis.NonBasis()[iv];
        if (inb < n) {
            for (SvIterator el(model_->GetCol(inb)); el; ++el) {
                res[el.index()] += el.value() * x[iv];
            }
        } else {
            res[inb - n] += x[iv];
        }
    }
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
