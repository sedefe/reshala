#pragma once

#include "reshala/linalg/dense_matrix.h"
#include "reshala/model/milp_model.h"
#include "reshala/model/solution.h"

namespace reshala {

struct DsState {
    std::vector<Index> basis;
    std::vector<Index> non_basis;
    std::vector<Index> index2nb;

    DenseVector c_n;
    DenseVector x_b;
    DenseVector x_n;
    std::vector<int8_t> d_p;

    DenseMatrix Binv;
};

class DualSimplex {
   public:
    DualSimplex(MilpModel& model);
    Solution Solve(bool warm);

    Index GetNIter() const { return n_iter; }

    DsState Store() const;
    void Restore(const DsState& state);

   private:
    MilpModel& model_;

    Index m, n;
    Index n_iter = 0;

    std::vector<Index> basis;
    std::vector<Index> non_basis;
    std::vector<Index> index2nb;

    DenseVector c_n;
    DenseVector x_b;
    DenseVector x_n;
    DenseVector e_p;
    DenseVector a_p;
    DenseVector a_q;
    std::vector<int8_t> d_p;  // типы небазисных переменных

    int8_t s_p;
    Scalar primal_infeasibility;
    Scalar theta_p, theta_d;
    Scalar a_pq;
    Scalar c_j_entering;

    Index iv_leaving, iv_entering;

    DenseMatrix Binv;

    void Init();
    void Chuzr();
    void Btran();
    void Price();
    void Chuzc();
    void Ftran();
    void Update();

    inline void SetXnValue(Index iv);
    Domain initial_domain;
    void ForceBounds();
    void UnforceBounds();

    void DebugPrint();
};

}  // namespace reshala
