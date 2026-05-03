#include "reshala/model/milp_model.h"

namespace reshala {

std::ostream& operator<<(std::ostream& os, const MilpModel& model) {
    os << model.obj << std::endl;

    os << "Subject to\n";
    Index m = model.GetNRows();
    Index n = model.GetNCols();
    for (Index i = 0; i < m; i++) {
        const SparseVector& lhs = model.Ar.getRows()[i];
        const Bound& rhs = model.rhs[i];
        if (rhs.le != -kInf) {
            os << lhs << " >= " << rhs.le << std::endl;
        }
        if (rhs.ri != kInf) {
            os << lhs << " <= " << rhs.ri << std::endl;
        }
    }

    os << "Bounds\n";
    for (Index i = 0; i < n; i++) {
        const Bound & bnd = model.bounds[i];
        if (bnd.le != -kInf) {
            os << "x[" << i << "] >= " << bnd.le << std::endl;
        }
        if (bnd.ri != kInf) {
            os << "x[" << i << "] <= " << bnd.ri << std::endl;
        }
    } 

    os << "Generals\n";
    for (Index i = 0; i < n; i++) {
        if (model.integrality[i]) {
            os << i << " ";
        }
    } 
    os << "\n";
    os << "End\n";

    return os;
}

}  // namespace reshala
