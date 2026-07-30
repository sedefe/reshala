#include "reshala/cuts/generators/cmir.h"

namespace reshala {

void CmirCg::Generate(const Solution& sol, std::vector<Cut>& dst) {
    Index m = model_.GetNCons();
    Index n = model_.GetNVars();
    const Index max_support = std::max(2, Index(kMaxRelSupport * n));
    x_scaled = ds_.GetX();

    for (Index ic = 0; ic < m; ic++) {
        if (!PrepareRow(ic)) continue;
        DoCut(dst);

        Cut cut(CutType::kCmir, lhs, rhs);
        // std::cout << "\tcut at row" << ic << ": " << cut;
        if (cut.IsViolated(sol.x) and lhs.Size() <= max_support) {
            dst.push_back(cut);
        }
    }
}

bool CmirCg::PrepareRow(Index ic) {
    Index m = model_.GetNCons();
    Index n = model_.GetNVars();
    const MilpModel& model_scaled = ds_.GetModel();

    Index ib = ds_.GetBasis().Basis()[ic];
    if (ib >= n) return false;                           // slack
    if (!model_scaled.GetIntegrality(ib)) return false;  // continuous

    // apply basic col scaling
    Scalar xb = x_scaled[ib];
    if (IsZero(MinFraction(xb))) return false;

    DenseVector lhs_dense;
    ds_.GetBasicRow(ic, lhs_dense);  // Btran+Price in scaled space

    lhs = SparseVector(lhs_dense);
    lhs.SetDim(n + m);
    for (MutableSvIterator el(lhs); el; ++el) {
        el.indexRef() = ds_.GetBasis().NonBasis()[el.index()];
    }
    lhs.Push(ib, 1.0);
    lhs.Sort();

    rhs = 0;

    return true;
}

void CmirCg::DoCut(std::vector<Cut>& dst) {
    Index m = model_.GetNCons();
    Index n = model_.GetNVars();
    const MilpModel& model_scaled = ds_.GetModel();
    std::vector<bool> sides(lhs.Size());

    // Displacement: x <- l+d or x <- u-d
    for (Index i = 0; i < lhs.Size(); i++) {
        Index iv = lhs.indices()[i];
        Scalar v = lhs.values()[i];
        const Bounds& bnd = model_scaled.GetBounds(iv);
        if (x_scaled[iv] - bnd.le > bnd.ri - x_scaled[iv]) {  // ri
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

        if (model_scaled.GetIntegrality(el.index())) {
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
        if (sides[i]) {  // ri
            lhs.values()[i] = -v;
            rhs -= v * model_scaled.GetBounds(iv).ri;
        } else {  // le
            rhs += v * model_scaled.GetBounds(iv).le;
        }
    }

    // Eliminate slacks. Todo: don't do this
    auto lhs_copy = lhs;
    {
        for (SvIterator el(lhs); el; ++el) {
            if (el.index() >= n) {
                lhs_copy = lhs_copy - el.value() * model_scaled.GetRow(el.index() - n);
                // lhs_copy.EraseIndex(el.index());  // Todo use EraseOffset()
            }
        }
        lhs_copy.SetDim(model_.GetNVars());
    }

    // Convert to our cuts format & unscale a_ij
    for (MutableSvIterator el(lhs_copy); el; ++el) {
        el.valueRef() = std::ldexp(el.value(), ds_.GetScaling().col[el.index()]);
    }

    std::swap(lhs, lhs_copy);
}

}  // namespace reshala
