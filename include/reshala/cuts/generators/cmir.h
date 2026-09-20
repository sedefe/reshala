#pragma once

#include "reshala/cuts/generators/abstract_cg.h"
#include "reshala/presolve/presolve.h"

namespace reshala {

class CmirCg : public AbstractCg {
   public:
    CmirCg(MilpModel& model, const DualSimplex& ds) : AbstractCg("Cmir", model, ds) {}

    void Generate(const Solution& sol, std::vector<Cut>& dst) override;

   private:
    DenseVector x;

    bool PrepareRow(Index ic, SparseVector& lhs);
    void DoCut(SparseVector& lhs, Scalar& rhs);
};

}  // namespace reshala
