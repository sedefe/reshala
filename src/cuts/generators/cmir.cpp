#include "reshala/cuts/generators/cmir.h"

namespace reshala {

void CmirCg::Generate(const Solution& sol, std::vector<Cut>& dst) {
    MilpModel& model = ctx_.GetModel();

    Index m = model.GetNCons();
    Index n = model.GetNVars();
    const Index max_support = std::max(2, Index(kMaxRelSupport * n));

    x = sol.x;
    x.resize(m + n);
    auto slacks = sol.slacks;
    std::copy(slacks.begin(), slacks.end(), x.begin() + n);

    SparseVector lhs(n);
    Scalar rhs;

    for (Index ic = 0; ic < m; ic++) {
        if (!PrepareRow(ic, lhs)) continue;
        DoCut(lhs, rhs);

        Cut cut(CutType::kCmir, lhs, rhs);
        // std::cout << "\tcut at row" << ic << ": " << cut;
        if (cut.IsViolated(sol.x) and lhs.Size() <= max_support) {
            dst.push_back(cut);
        }
    }
}

bool CmirCg::PrepareRow(Index ic, SparseVector& lhs) {
    MilpModel& model = ctx_.GetModel();
    DualSimplex& ds = ctx_.GetDs();

    Index m = model.GetNCons();
    Index n = model.GetNVars();

    Index ib = ds.GetBasis().Basis()[ic];
    if (ib >= n) return false;                    // slack
    if (!model.GetIntegrality(ib)) return false;  // continuous

    // apply basic col scaling
    Scalar xb = x[ib];
    if (IsZero(MinFraction(xb))) return false;

    DenseVector lhs_dense;
    ds.GetBasicRow(ic, lhs_dense);  // Btran+Price in scaled space

    lhs = SparseVector(lhs_dense);
    for (MutableSvIterator el(lhs); el; ++el) {
        el.indexRef() = ds.GetBasis().NonBasis()[el.index()];
    }

    // Unscale but keep ib's coeff as 1
    Scalar c = ds.GetScaling().col[ib];
    Index scale;
    for (MutableSvIterator el(lhs); el; ++el) {
        if (el.index() < n) {
            scale = ds.GetScaling().col[el.index()] - c;
        } else {
            scale = -ds.GetScaling().row[el.index() - n] - c;
        }
        el.valueRef() = std::ldexp(el.value(), scale);
    }
    lhs.Push(ib, 1.0);
    lhs.Sort();

    return true;
}

void CmirCg::DoCut(SparseVector& lhs, Scalar& rhs) {
    MilpModel& model = ctx_.GetModel();

    Index m = model.GetNCons();
    Index n = model.GetNVars();
    std::vector<bool> sides(lhs.Size());

    rhs = 0;
    // Displacement: x <- l+d or x <- u-d
    for (Index i = 0; i < lhs.Size(); i++) {
        Index iv = lhs.indices()[i];
        Scalar v = lhs.values()[i];

        Bounds bnd = (iv < n) ? model.GetBounds(iv)
                              : Bounds{-model.GetRhs(iv - n).ri, -model.GetRhs(iv - n).le};

        if (x[iv] - bnd.le > bnd.ri - x[iv]) {  // ri
            rhs -= v * bnd.ri;
            lhs.values()[i] = -v;
            sides[i] = true;
        } else {  // le
            rhs -= v * bnd.le;
            sides[i] = false;
        }
    }

    Scalar f = Fraction(rhs);

    // Generate cut coeffs
    // sum(aj dj) = b  =>  sum(alphaj dj) >= ceil(b) * frac(b)
    for (MutableSvIterator el(lhs); el; ++el) {
        Scalar r = Fraction(el.value());
        Scalar v = el.value();

        if (el.index() < n and model.GetIntegrality(el.index())) {
            if (r > f) {
                el.valueRef() = f * Ceil(el.value());
            } else {
                el.valueRef() = f * Ceil(el.value()) + r;
            }
        } else {
            if (el.value() < 0) {
                el.valueRef() = 0.0;
            } else {
                el.valueRef() = el.value();
            }
        }
    }
    rhs = Ceil(rhs) * f;

    // Backward substitution
    for (Index i = 0; i < lhs.Size(); i++) {
        Index iv = lhs.indices()[i];
        Scalar v = lhs.values()[i];

        Bounds bnd = (iv < n) ? model.GetBounds(iv)
                              : Bounds{-model.GetRhs(iv - n).ri, -model.GetRhs(iv - n).le};

        if (sides[i]) {  // ri
            lhs.values()[i] = -v;
            rhs -= v * bnd.ri;
        } else {  // le
            rhs += v * bnd.le;
        }
    }

    // Eliminate slacks. Todo: don't do this
    auto lhs_copy = lhs;
    {
        for (SvIterator el(lhs); el; ++el) {
            if (el.index() >= n) {
                lhs_copy = axpy(-el.value(), model.GetRow(el.index() - n), lhs_copy);
                lhs_copy.EraseIndex(el.index());  // Todo use EraseOffset()
            }
        }
        lhs_copy.SetDim(model.GetNVars());
    }

    std::swap(lhs, lhs_copy);
}

}  // namespace reshala
