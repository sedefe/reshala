#pragma once

#include "reshala/lina/lina.h"
#include "reshala/lp/scaling.h"
#include "reshala/model/milp_model.h"
#include "reshala/model/solution.h"

namespace reshala {

struct DsState {
    DenseVector c_n;
    DenseVector x_b;
    std::vector<int8_t> d_n;

    LpBasis basis;
    Lina lina;
};

struct DsStats {
    Index n_iter = 0;
};
std::ostream& operator<<(std::ostream& os, const DsStats& stats);

class DualSimplex {
    const Scalar kPivotTolerance = 1e-6;
    const Index kMaxLinaAge = 50;

   public:
    DualSimplex() {}
    void SetModel(MilpModel& model);
    Solution Solve(bool warm);

    inline const DsStats& GetStats() const { return stats; }
    inline const LpBasis& GetBasis() const { return basis; }
    inline const Lina& GetLina() const { return lina; }
    inline const Scaling& GetScaling() const { return scaling; }
    inline const std::vector<int8_t>& GetDn() const { return d_n; }

    DsState Store() const;
    void Restore(const DsState& state);

    void SetDomain(const Domain& domain) {
        model_orig_->SetDomain(domain);
        for (Index iv = 0; iv < domain.Size(); iv++) {
            SetBounds(iv, domain.GetBounds(iv));
        }
    }
    inline void SetBounds(Index iv, const Bounds& bnd) {
        model_orig_->SetBounds(iv, bnd);
        model_.SetBounds(
            iv, {std::ldexp(bnd.le, scaling.col[iv]), std::ldexp(bnd.ri, scaling.col[iv])});
    }

    void GetBasicRow(Index ic, DenseVector& res) const;

   private:
    MilpModel* model_orig_;
    MilpModel model_;  // Scaled
    Scaling scaling;

    DsStats stats;

    Index m, n;
    LpStatus status;
    LpBasis basis;

    DenseVector c_n;
    DenseVector x_b;
    DenseVector a_p;
    DenseVector a_q;
    std::vector<int8_t> d_n;  // типы небазисных переменных

    int8_t s_p;
    Scalar primal_infeasibility;
    Scalar theta_p, theta_d;
    Scalar a_pq;
    Scalar c_q;

    Index iv_leaving, iv_entering;

    Lina lina;

    void Init();
    void Chuzr();
    void Chuzc();
    void Ftran();
    void Update();
    void RebuildAll();

    Solution PrepareSolution() const;

    Scalar GetXnValue(Index iv) const;
    void MulNLeft(const DenseVector& x, DenseVector& res) const;
    void MulNRight(const DenseVector& x, DenseVector& res) const;

    void DebugPrint();
};

}  // namespace reshala
