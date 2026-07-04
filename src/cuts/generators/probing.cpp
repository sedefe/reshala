#include "reshala/cuts/generators/probing.h"

namespace reshala {

Index ProbingCg::Generate(const Solution& sol, std::vector<Cut>& dst) {
    Index dim = model_.GetNVars();
    for (const auto& impl : impls_) {
        Index iv1 = impl.x_ind;
        Index iv2 = impl.y_ind;

        bool x1 = impl.x_value;

        Scalar b = impl.b;
        bool left = impl.left;

        const Bounds& bnd = model_.GetBounds(iv2);
        if (left) {
            if (x1 == false) {
                // x = 0 => y >= b
                // x = 1 => y >= bnd.le
                // y >= b + (bnd.le-b) * x
                // (b - bnd.le) x + y >= b
                Cut cut(SparseVector(dim, {iv1, iv2}, {b - bnd.le, 1.0}, false), b);
                auto v = cut.GetViol(sol.x);
                if (StrongGt(v, 0)) {
                    dst.push_back(cut);
                }
            } else {
                // x = 0 => y >= bnd.le
                // x = 1 => y >= b
                // y >= bnd.le + (b - bnd.le) * x
                // (bnd.le - b) x + y >= bnd.le
                Cut cut(SparseVector(dim, {iv1, iv2}, {bnd.le - b, 1.0}, false), bnd.le);
                auto v = cut.GetViol(sol.x);
                if (StrongGt(v, 0)) {
                    dst.push_back(cut);
                }
            }
        } else {
            if (x1 == false) {
                // x = 0 => y <= b
                // x = 1 => y <= bnd.ri
                // y <= b + (bnd.ri - b) * x
                // (bnd.ri - b) x - y >= - b
                Cut cut(SparseVector(dim, {iv1, iv2}, {bnd.ri - b, -1.0}, false), -b);
                auto v = cut.GetViol(sol.x);
                if (StrongGt(v, 0)) {
                    dst.push_back(cut);
                }
            } else {
                // x = 0 => y <= bnd.ri
                // x = 1 => y <= b
                // y <= bnd,ri + (b - bnd.ri) * x
                // (b - bnd.ri) x - y >= - bnd.ri
                Cut cut(SparseVector(dim, {iv1, iv2}, {b - bnd.ri, -1.0}, false), -bnd.ri);
                auto v = cut.GetViol(sol.x);
                if (StrongGt(v, 0)) {
                    dst.push_back(cut);
                }
            }
        }
    }

    return 0;  // Todo
}

}  // namespace reshala
