#pragma once

#include "reshala/lp/dual_simplex.h"
#include "reshala/model/milp_model.h"
#include "reshala/utils.h"

namespace reshala {

struct Node {
    Node() {}
    Node(Index l, const Solution& s, const Domain& d, const DsState& st)
        : level(l), sol(s), domain(d), ds_state(st) {}
    Index level;
    Solution sol;
    Domain domain;

    DsState ds_state;
};

class MipState {
   public:
    MipState(const MilpModel& model, const Solution best_sol = InfeasibleSolution(),
             Scalar dual = -kInf)
        : model_(model), best_sol_(best_sol), dual_(dual) {
        int_obj_ = ObjIsInteger();
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
    bool int_obj_;
    Scalar gap_;
    Solution best_sol_;

    bool ObjIsInteger() {
        for (Index iv = 0; iv < model_.GetNVars(); ++iv) {
            if (model_.GetIntegrality(iv)) {
                if (GetFraction(model_.GetObj().coefficients[iv]) != 0.0) {
                    return false;
                }
            } else if (model_.GetObj().coefficients[iv] != 0.0) {
                return false;
            }
        }
        return true;
    }

    void Recalc() {
        Scalar primal = best_sol_.y;
        Scalar min_diff = std::max(std::abs(primal) * kMipGap, int_obj_ ? 1.0 - kEpsZero : 0.0);
        cutoff_ = best_sol_.status == LpStatus::kInfeasible ? kInf : primal - min_diff;
        gap_ = IsZero(primal) ? (IsZero(dual_) ? 0 : 1)
                              : (primal == kInf ? kInf : (primal - dual_) / std::abs(primal));
    }
};

}  // namespace reshala
