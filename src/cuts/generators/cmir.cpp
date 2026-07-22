#include "reshala/cuts/generators/cmir.h"

namespace reshala {

void CmirCg::Generate(const Solution& sol, std::vector<Cut>& dst) {
    Index m = model_.GetNCons();
    Index n = model_.GetNVars();
    const MilpModel& model_scaled = ds_.GetModel();
    const Index max_support = std::max(2, Index(kMaxRelSupport * n));
    const auto& d_n = ds_.GetDn();

    for (Index ic = 0; ic < m; ic++) {
        Index ib = ds_.GetBasis().Basis()[ic];
        if (ib >= n) continue;  // slack
        if (!model_scaled.GetIntegrality(ib)) continue;

        // apply basic col scaling
        Scalar xb = std::ldexp(sol.x[ib], ds_.GetScaling().col[ib]);
        if (IsZero(MinFraction(xb))) continue;

        DenseVector x;
        ds_.GetBasicRow(ic, x);  // Btran+Price in scaled space

        SparseVector lhs(x);
        lhs.SetDim(model_scaled.GetNVars());
        if (lhs.Size() > max_support) continue;

        Scalar f = Fraction(xb);

        // Displacement: x <- l+d or x <- u-d
        for (MutableSvIterator el(lhs); el; ++el) {
            if (d_n[el.index()] < 0) {
                el.valueRef() = -el.value();
            }
        }

        // Generate cut coeffs
        // xb + sum(aj dj) = b  =>  xb + sum(alphaj dj) <= floor(b)
        for (MutableSvIterator el(lhs); el; ++el) {
            Index inb = ds_.GetBasis().NonBasis()[el.index()];
            Scalar r = Fraction(el.value());

            if (inb < n and model_scaled.GetIntegrality(inb)) {
                el.valueRef() = Floor(el.value()) + std::max(0.0, r - f) / (1 - f);
            } else {
                if (el.value() > 0) {
                    el.valueRef() = 0.0;
                } else {
                    el.valueRef() = el.value() / (1 - f);
                }
            }
        }

        // Backward substitution
        Scalar b = Floor(xb);
        for (MutableSvIterator el(lhs); el; ++el) {
            Index inb = ds_.GetBasis().NonBasis()[el.index()];

            if (d_n[el.index()] < 0) {
                el.valueRef() = -el.value();
                b += el.value() * model_scaled.GetBounds(inb).ri;
            } else {
                b += el.value() * model_scaled.GetBounds(inb).le;
            }
            el.indexRef() = inb;  // replace nonbasic indices with original ones
        }

        lhs.Push(ib, 1.0);  // Add x_b
        lhs.Sort();

        auto lhs_copy = lhs;
        {  // Eliminate slacks. Todo: don't do this
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
            el.valueRef() = -std::ldexp(el.value(), ds_.GetScaling().col[el.index()]);
        }

        Cut cut(CutType::kCmir, lhs_copy, -b);
        // std::cout << "\tcut at row" << ic << ": " << cut;
        if (cut.IsViolated(sol.x)) {
            dst.push_back(cut);
        }
    }
}

}  // namespace reshala
