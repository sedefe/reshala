#include "reshala/lp/scaling.h"

namespace reshala {

std::ostream& operator<<(std::ostream& os, const ScaleReport& rep) {
    os << "FrobNorm: " << rep.frob_norm << ", MaxMinRatio: " << rep.max_min_ratio;
    return os;
}

std::ostream& operator<<(std::ostream& os, const ScaleStats& stats) {
    os << "Scaling stats:\n"
       << "\tBefore: " << stats.before << "\n"
       << "\tAfter:  " << stats.after << "\n";
    return os;
}

ScaleReport::ScaleReport(const SparseRowMatrix& A) {
    frob_norm = 0.0;
    Index nnz = A.GetNnz();
    Scalar mx = -kInf, mn = kInf, val;
    for (Index ic = 0; ic < A.GetNRows(); ic++) {
        if (A.GetRow(ic).Size() == 0) continue;
        for (SvIterator el(A.GetRow(ic)); el; ++el) {
            val = std::abs(el.value());
            frob_norm += val * val;
            mx = std::max(mx, val);
            mn = std::min(mn, val);
        }
    }
    max_min_ratio = mx / mn;
}

void Scaling::ScaleModel(MilpModel& model) {
    Index m = model.GetNCons();
    Index n = model.GetNVars();

    row.assign(m, 0.0);
    col.assign(n, 0.0);

    stats.before = ScaleReport(model.GetAr());

    // Init scaling
    Index exp;
    Index sum;
    Index kMaxIters = 5;

    for (Index i = 0; i < kMaxIters; i++) {
        for (Index ic = 0; ic < m; ic++) {
            if (model.GetRow(ic).Size() == 0) continue;
            sum = 0;
            for (SvIterator el(model.GetRow(ic)); el; ++el) {
                std::frexp(el.value(), &exp);
                sum += exp - col[el.index()];
            }
            row[ic] = Scalar(sum) / model.GetRow(ic).Size() + 0.5;
        }

        for (Index iv = 0; iv < n; iv++) {
            if (model.GetAc().GetCol(iv).Size() == 0) continue;
            sum = 0;
            for (SvIterator el(model.GetCol(iv)); el; ++el) {
                std::frexp(el.value(), &exp);
                sum += exp - row[el.index()];
            }
            col[iv] = Scalar(sum) / model.GetCol(iv).Size() + 0.5;
        }
    }

    // Apply
    for (Index iv = 0; iv < n; iv++) {
        for (MutableSvIterator el(model.GetCol(iv)); el; ++el) {
            el.valueRef() = std::ldexp(el.value(), -row[el.index()] - col[iv]);
        }
        model.GetObj().coefficients[iv] = std::ldexp(model.GetObj().coefficients[iv], -col[iv]);
        const Bounds& bnd = model.GetBounds(iv);
        model.SetBounds(iv, {std::ldexp(bnd.le, col[iv]), std::ldexp(bnd.ri, col[iv])});
    }
    for (Index ic = 0; ic < m; ic++) {
        for (MutableSvIterator el(model.GetRow(ic)); el; ++el) {
            el.valueRef() = std::ldexp(el.value(), -row[ic] - col[el.index()]);
        }
        Bounds& rhs = model.GetRhs(ic);
        rhs = {std::ldexp(rhs.le, -row[ic]), std::ldexp(rhs.ri, -row[ic])};
    }

    stats.after = ScaleReport(model.GetAr());
}

}  // namespace reshala
