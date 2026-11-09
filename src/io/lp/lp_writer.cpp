#include "reshala/io/lp/lp_writer.h"

#include "reshala/io/io.h"
#include "reshala/presolve/presolve.h"

namespace reshala {

LpWriter::LpWriter(Reshala& ctx) : ctx_(ctx) {}

void LpWriter::PrintObj(std::ostream& os, const Objective& obj) const {
    const auto& names = ctx_.GetIo().GetNames().vars;
    const auto& orig_var_idx = ctx_.GetPresolver().GetTracker().GetOrigVarIdx();

    Scalar sense = 0;
    bool first_coeff = true;
    if (obj.orig_sense == Sense::kMin) {
        os << "Minimize\n";
        sense = 1.;
    } else {
        os << "Maximize\n";
        sense = -1.;
    }
    for (Index i = 0; i < obj.coefficients.size(); i++) {
        if (obj.coefficients[i] == 0) continue;
        const auto& name = names.GetName(orig_var_idx[i]);
        Scalar v = sense * obj.coefficients[i] * obj.mult;

        if (first_coeff) {
            os << (v > 0 ? "" : "- ");
        } else {
            os << (v > 0 ? " + " : " - ");
        }
        first_coeff = false;

        v = std::abs(v);
        if (v == 1) {
            os << " " << name;
        } else {
            os << v << " " << name;
        }
    }
    if (obj.c0 != 0) {
        if (first_coeff) {
            os << (obj.c0 > 0 ? "" : "- ");
        } else {
            os << (obj.c0 > 0 ? " + " : " - ");
        }
        os << std::abs(obj.c0 * obj.mult);
    }

    os << "\n";
}

void LpWriter::PrintCon(std::ostream& os, const SparseVector& lhs, const std::string& exp,
                        Scalar rhs) const {
    const auto& names = ctx_.GetIo().GetNames().vars;
    const auto& orig_var_idx = ctx_.GetPresolver().GetTracker().GetOrigVarIdx();

    for (Index i = 0; i < lhs.Size(); ++i) {
        const auto& name = names.GetName(orig_var_idx[lhs.indices()[i]]);
        Scalar v = lhs.values()[i];

        if (i == 0) {
            os << (v >= 0 ? "" : "- ");
        } else {
            os << (v >= 0 ? " + " : " - ");
        }
        v = std::abs(v);
        if (v == 1) {
            os << name;
        } else {
            os << v << " " << name;
        }
    }

    os << " " << exp << " " << rhs << "\n";
}

void LpWriter::PrintModel(std::ostream& os, const MilpModel& model) const {
    const auto& names = ctx_.GetIo().GetNames().vars;
    const auto& orig_var_idx = ctx_.GetPresolver().GetTracker().GetOrigVarIdx();

    PrintObj(os, model.GetObj());

    os << "Subject to\n";
    Index m = model.GetNCons();
    Index n = model.GetNVars();
    for (Index i = 0; i < m; i++) {
        const auto& lhs = model.GetRow(i);
        const auto& rhs = model.GetRhs(i);

        if (rhs.le == rhs.ri) {
            PrintCon(os, lhs, "=", rhs.le);
        } else {
            if (rhs.le != -kInf) PrintCon(os, lhs, ">=", rhs.le);
            if (rhs.ri != kInf) PrintCon(os, lhs, "<=", rhs.ri);
        }
    }

    std::vector<Index> binaries, generals;
    os << "Bounds\n";
    for (Index iv = 0; iv < n; iv++) {
        const auto& name = names.GetName(orig_var_idx[iv]);
        const auto& bnd = model.GetBounds(iv);

        if (model.GetType(iv) == BndType::kFixed) {
            os << name << " = " << (bnd.le + bnd.ri) / 2 << std::endl;
        } else {
            if (model.GetIntegrality(iv)) {
                if (bnd.le == 0 and bnd.ri == 1) {
                    binaries.push_back(iv);
                    continue;
                } else {
                    generals.push_back(iv);
                }
            }
            if (bnd.le != 0.0) os << name << " >= " << bnd.le << std::endl;
            if (bnd.ri != kInf) os << name << " <= " << bnd.ri << std::endl;
        }
    }

    if (binaries.size()) {
        os << "Binaries\n";
        for (const auto& iv : binaries) os << names.GetName(orig_var_idx[iv]) << " ";
        os << "\n";
    }
    if (generals.size()) {
        os << "Generals\n";
        for (const auto& iv : generals) os << names.GetName(orig_var_idx[iv]) << " ";
        os << "\n";
    }

    os << "End\n";
}

}  // namespace reshala
