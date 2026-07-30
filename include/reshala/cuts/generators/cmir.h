#pragma once

#include "reshala/cuts/generators/abstract_cg.h"
#include "reshala/presolve/presolve.h"

namespace reshala {

class CmirCg : public AbstractCg {
   public:
    CmirCg(MilpModel& model, DualSimplex& ds)
        : AbstractCg("Cmir", model, ds), lhs(model.GetNVars()) {}

    void Generate(const Solution& sol, std::vector<Cut>& dst) override;

   private:
    bool PrepareRow(Index ic);
    void DoCut(std::vector<Cut>& dst);
    SparseVector lhs;
    Scalar rhs;
    DenseVector x;
};

}  // namespace reshala
