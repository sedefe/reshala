#pragma once

#include <numeric>

#include "reshala/lp/dual_simplex.h"
#include "reshala/model/milp_model.h"
#include "reshala/utils.h"

namespace reshala {

struct Node {
    Node() {}
    Node(Index l, const Solution& s, const Domain& d, const DsState& st)
        : level(l), sol(s), domain(d), ds_state(st) {
    }
    Index level;
    Solution sol;
    Domain domain;

    DsState ds_state;
};

class MipState {
   public:
    MipState(const MilpModel& model, const Solution best_sol = {LpStatus::kInfeasible, kInf, {}},
             Scalar dual = -kInf)
        : model_(model), best_sol_(best_sol), dual_(dual) {
        obj_gcd_ = GetObjectiveGcd();
        Recalc();
    }
    const Solution& GetBestSol() const { return best_sol_; }
    Scalar GetPrimal() const { return best_sol_.y; }
    Scalar GetDual() const { return dual_; }
    Scalar GetGap() const { return gap_; }
    Scalar GetCutoff() const { return cutoff_; };

    bool TestPrimal(const Solution& sol) {
        if (sol.status != LpStatus::kOptimal) return false;
        if (best_sol_.y <= sol.y) return false;
        if (!model_.IsIntegerFeasible(sol.x)) return false;
        best_sol_ = sol;
        Recalc();
        return true;
    }
    void UpdDual(Scalar dual) {
        dual_ = dual;
        Recalc();
    }
    bool Converged() const { return cutoff_ <= dual_; }

   private:
    const MilpModel& model_;
    Scalar dual_;
    Scalar cutoff_;
    Scalar obj_gcd_;
    Scalar gap_;
    Solution best_sol_;

    Index GetObjectiveGcd() {
        Index res = 0;
        for (Index iv = 0; iv < model_.GetNVars(); ++iv) {
            if (model_.GetIntegrality(iv)) {
                if (GetFraction(model_.GetObj().coefficients[iv]) == 0.0) {
                    res = std::gcd(res, Index(model_.GetObj().coefficients[iv]));
                } else {
                    return 0;
                }
            } else if (model_.GetObj().coefficients[iv] != 0.0) {
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
