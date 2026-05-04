#pragma once

#include "reshala/linalg/sparse_matrix.h"
#include "reshala/model/domain.h"
#include "reshala/model/objective.h"

namespace reshala {

class MilpModel {
   public:
    MilpModel() {}

    Index GetNRows() const { return Ac.getNCols(); }
    Index GetNCols() const { return Ac.getNRows(); }

    const Objective& GetObj() const { return obj; }
    Objective& GetObj() { return obj; }

    const SparseColMatrix& GetAc() const { return Ac; }
    SparseColMatrix& GetAc() { return Ac; }

    const SparseRowMatrix& GetAr() const { return Ar; }
    SparseRowMatrix& GetAr() { return Ar; }

    const std::vector<Bounds>& GetRhs() const { return rhs; }
    std::vector<Bounds>& GetRhs() { return rhs; }

    const Domain& GetVars() const { return vars; }
    Domain& GetVars() { return vars; }

    void Resize(Index m, Index n) {
        obj.coefficients.resize(n);
        Ac.Resize(m, n);
        Ar.Resize(m, n);
        rhs.resize(m);
        vars.Resize(n);
    }

    void PrepareConstraint(const SparseVector& lhs, const Bounds& rhs) {
        this->Ar.GetRows().push_back(lhs);
        this->rhs.push_back(rhs);
    }

    friend std::ostream& operator<<(std::ostream& os, const MilpModel& model);

   private:
    Objective obj;
    SparseColMatrix Ac;
    SparseRowMatrix Ar;
    std::vector<Bounds> rhs;
    Domain vars;
};

}  // namespace reshala
