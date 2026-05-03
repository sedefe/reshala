#pragma once

#include "reshala/linalg/dense_matrix.h"
#include "reshala/model/milp_model.h"
#include "reshala/model/solution.h"

namespace reshala {

class DualRevisedSimplex {
   public:
    DualRevisedSimplex(MilpModel& model_);
    Solution solve();

   private:
    MilpModel& model;

    Index m, n;
    Index n_iter = 0;

    std::vector<Index> basis;
    std::vector<Index> non_basis;

    DenseVector c_b;
    DenseVector c_n;
    DenseVector x_b;
    DenseVector x_n;
    DenseVector e_p;
    DenseVector a_p;
    DenseVector a_q;

    Index iv_leaving, iv_entering;

    DenseMatrix Binv;

    void init();
    void chuzr();
    void btran();
    void price();
    void chuzc();
    void ftran();
    void update();
};

}  // namespace reshala
