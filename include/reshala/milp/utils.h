#pragma once

#include <numeric>

#include "reshala/model/milp_model.h"
#include "reshala/utils.h"

namespace reshala {

struct Node {
    Node(const Solution& s, const Domain& d) : sol(s), domain(d) {}
    Solution sol;
    Domain domain;
};

class MipState {
   public:
    MipState(Scalar dual, Scalar primal, const MilpModel& model) {
        primal_ = primal;
        dual_ = dual;
        obj_gcd_ = GetObjectiveGcd(model);
        Recalc();
    }
    Scalar GetPrimal() const { return primal_; }
    Scalar GetDual() const { return dual_; }
    Scalar GetGap() const { return gap_; }
    Scalar GetCutoff() const { return cutoff_; };

    void UpdPrimal(Scalar primal) {
        primal_ = primal;
        Recalc();
    }
    void UpdDual(Scalar dual) {
        dual_ = dual;
        Recalc();
    }
    bool Converged() const { return cutoff_ <= dual_; }

   private:
    Scalar dual_;
    Scalar cutoff_;
    Scalar obj_gcd_;
    Scalar primal_;
    Scalar gap_;

    Index GetObjectiveGcd(const MilpModel& model) {
        Index res = 0;
        for (Index iv = 0; iv < Index(model.GetNVars()); ++iv) {
            if (model.GetIntegrality(iv)) {
                if (GetFraction(model.GetObj().coefficients[iv]) == 0.0) {
                    res = std::gcd(res, Index(model.GetObj().coefficients[iv]));
                } else {
                    return 0;
                }
            } else if (model.GetObj().coefficients[iv] != 0.0) {
                return 0;
            }
        }
        return res;
    }

    void Recalc() {
        cutoff_ = primal_ == kInf ? kInf
                                  : primal_ - std::max(abs(primal_) * kMipGap, obj_gcd_ - kEpsZero);
        gap_ = primal_ == 0 ? dual_ == 0 ? 0 : kInf : abs(1 - dual_ / primal_);
    }
};

}  // namespace reshala
