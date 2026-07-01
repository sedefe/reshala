#pragma once

#include "reshala/lina/lina.h"
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
    Scalar kPivotTolerance = 1e-6;

   public:
    DualSimplex() {}
    void SetModel(MilpModel& model);
    Solution Solve(bool warm);

    inline const DsStats& GetStats() const { return stats; }

    DsState Store() const;
    void Restore(const DsState& state);

   private:
    MilpModel* model_;

    DsStats stats;

    Index m, n;

    LpBasis basis;

    DenseVector c_n;
    DenseVector x_b;
    DenseVector e_p;
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
    void Btran();
    void Price();
    void Chuzc();
    void Ftran();
    void Update();

    inline Scalar GetXnValue(Index iv);
    Domain initial_domain;
    void ForceBounds();
    void UnforceBounds();

    void DebugPrint();
};

}  // namespace reshala
