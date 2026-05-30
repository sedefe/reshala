#pragma once

#include "reshala/model/milp_model.h"
#include "reshala/utils.h"

namespace reshala {

class Activity {
    Bounds finite;    // sum of finite parts (never +/-inf)
    Index n_neg_inf;  // count of terms with -inf lower bound
    Index n_pos_inf;  // count of terms with +inf upper bound

   public:
    Activity() : finite{0, 0}, n_neg_inf(0), n_pos_inf(0) {}

    void AddTerm(Scalar coef, const Bounds& bnd) {
        if (StrongGt(coef, 0)) {
            finite.le += coef * (bnd.le == -kInf ? 0 : bnd.le);
            finite.ri += coef * (bnd.ri == kInf ? 0 : bnd.ri);
            n_neg_inf += (bnd.le == -kInf);
            n_pos_inf += (bnd.ri == kInf);
        } else if (StrongLt(coef, 0)) {
            finite.le += coef * (bnd.ri == kInf ? 0 : bnd.ri);
            finite.ri += coef * (bnd.le == -kInf ? 0 : bnd.le);
            n_neg_inf += (bnd.ri == kInf);
            n_pos_inf += (bnd.le == -kInf);
        }
    }

    void RmTerm(Scalar coef, const Bounds& bnd) {
        if (StrongGt(coef, 0)) {
            finite.le -= coef * (bnd.le == -kInf ? 0 : bnd.le);
            finite.ri -= coef * (bnd.ri == kInf ? 0 : bnd.ri);
            n_neg_inf -= (bnd.le == -kInf);
            n_pos_inf -= (bnd.ri == kInf);
        } else if (StrongLt(coef, 0)) {
            finite.le -= coef * (bnd.ri == kInf ? 0 : bnd.ri);
            finite.ri -= coef * (bnd.le == -kInf ? 0 : bnd.le);
            n_neg_inf -= (bnd.ri == kInf);
            n_pos_inf -= (bnd.le == -kInf);
        }
    }

    Bounds GetRange() const {
        Bounds result;
        result.le = (n_neg_inf > 0) ? -kInf : finite.le;
        result.ri = (n_pos_inf > 0) ? kInf : finite.ri;
        return result;
    }

    void Scale(Scalar x) {
        finite.le *= x;
        finite.ri *= x;
    }
};

}  // namespace reshala
