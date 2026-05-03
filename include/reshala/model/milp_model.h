#pragma once

#include "reshala/linalg/sparse_matrix.h"
#include "reshala/model/objective.h"

namespace reshala {

struct Bound {
    Scalar le = -kInf;
    Scalar ri = kInf;
};

class MilpModel {
   public:
    MilpModel() {}

    Index GetNRows() const { return Ac.GetNCols(); }
    Index GetNCols() const { return Ac.GetNRows(); }

    const Objective& getObj() const { return obj; }
    Objective& getObj() { return obj; }

    const SparseColMatrix& getAc() const { return Ac; }
    SparseColMatrix& getAc() { return Ac; }

    const SparseRowMatrix& getAr() const { return Ar; }
    SparseRowMatrix& getAr() { return Ar; }

    const std::vector<Bound>& getRhs() const { return rhs; }
    std::vector<Bound>& getRhs() { return rhs; }

    const std::vector<Bound>& getBounds() const { return bounds; }
    std::vector<Bound>& getBounds() { return bounds; }

    const std::vector<bool>& getIntegrality() const { return integrality; }
    std::vector<bool>& getIntegrality() { return integrality; }

    void resize(Index m, Index n) {
        obj.coefficients.resize(n);
        Ac.resize(m, n);
        Ar.resize(m, n);
        rhs.resize(m);
        bounds.resize(n);
        integrality.resize(n);
    }

    void PrepareConstraint(const SparseVector& lhs, const Bound& rhs) {
        this->Ar.getRows().push_back(lhs);
        this->rhs.push_back(rhs);
    }

    friend std::ostream& operator<<(std::ostream& os, const MilpModel& model);

   private:
    Objective obj;
    SparseColMatrix Ac;
    SparseRowMatrix Ar;
    std::vector<Bound> rhs;
    std::vector<Bound> bounds;
    std::vector<bool> integrality;
};

}  // namespace reshala
