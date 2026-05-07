#pragma once

#include "reshala/linalg/dense_matrix.h"
#include "reshala/model/milp_model.h"
#include "reshala/model/solution.h"

namespace reshala {

class DualSimplex {
   public:
    DualSimplex(MilpModel& model);
    Solution Solve();

   private:
    MilpModel& model_;

    Index m, n;
    Index n_iter = 0;

    std::vector<Index> basis;
    std::vector<Index> non_basis;
    std::vector<Index> index2nb;

    DenseVector c_b;
    DenseVector c_n;
    DenseVector x_b;
    DenseVector x_n;
    DenseVector e_p;
    DenseVector a_p;
    DenseVector a_q;

    int8_t s_p;
    Scalar primal_infeasibility;
    Scalar theta_p;
    Scalar theta_d;

    Index iv_leaving, iv_entering;

    DenseMatrix Binv;

    void Init();
    void Chuzr();
    void Btran();
    void Price();
    void Chuzc();
    void Ftran();
    void Update();

    Domain initial_domain;
    void ForceBounds();
    void UnforceBounds();

    void DebugPrint();
};

}  // namespace reshala
