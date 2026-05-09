#pragma once

#include <numeric>

#include "reshala/model/milp_model.h"
#include "reshala/utils.h"

namespace reshala {

struct Node {
    Node() {}
    Node(Index l, const Solution& s, const Domain& d) : level(l), sol(s), domain(d) {}
    Index level;
    Solution sol;
    Domain domain;
};

class MipState {
   public:
    MipState(const MilpModel& model) {
        best_sol_ = {LpStatus::kInfeasible, kInf, {}};
        dual_ = -kInf;
        obj_gcd_ = GetObjectiveGcd(model);
        Recalc();
    }
    const Solution& GetBestSol() const { return best_sol_; }
    Scalar GetPrimal() const { return best_sol_.y; }
    Scalar GetDual() const { return dual_; }
    Scalar GetGap() const { return gap_; }
    Scalar GetCutoff() const { return cutoff_; };

    bool TestPrimal(const Solution& sol) {
        if (best_sol_.y > sol.y) {
            best_sol_ = sol;
            Recalc();
            return true;
        }
        return false;
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
    Scalar gap_;
    Solution best_sol_;

    Index GetObjectiveGcd(const MilpModel& model) {
        Index res = 0;
        for (Index iv = 0; iv < model.GetNVars(); ++iv) {
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
        Scalar primal = best_sol_.y;
        cutoff_ = best_sol_.status == LpStatus::kInfeasible
                      ? kInf
                      : primal - std::max(abs(primal) * kMipGap, obj_gcd_ - kEpsZero);
        gap_ = primal == 0 ? (dual_ == 0 ? 0 : kInf) : std::abs(1 - dual_ / primal);
    }
};

}  // namespace reshala
