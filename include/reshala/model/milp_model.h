#pragma once

#include "reshala/linalg/sparse_matrix.h"
#include "reshala/model/domain.h"
#include "reshala/model/objective.h"
#include "reshala/model/solution.h"

namespace reshala {

class MilpModel {
   public:
    MilpModel() {}

    Index GetNCons() const { return Ac_.getNRows(); }
    Index GetNVars() const { return Ac_.getNCols(); }

    const Objective& GetObj() const { return obj_; }
    Objective& GetObj() { return obj_; }

    const SparseColMatrix& GetAc() const { return Ac_; }
    SparseColMatrix& GetAc() { return Ac_; }

    const SparseRowMatrix& GetAr() const { return Ar_; }
    SparseRowMatrix& GetAr() { return Ar_; }

    const std::vector<Bounds>& GetRhs() const { return rhs_; }
    std::vector<Bounds>& GetRhs() { return rhs_; }

    const Domain& GetVars() const { return domain_; }
    Domain& GetVars() { return domain_; }

    inline const Domain GetDomain() const { return domain_; }
    inline void SetDomain(const Domain& domain) { domain_ = domain; }
    inline const Bounds& GetBounds(Index iv) const { return domain_.GetBounds(iv); }
    inline void SetBounds(Index iv, const Bounds& bnd) { domain_.SetBounds(iv, bnd); }
    inline const VarType& GetType(Index iv) const { return domain_.GetType(iv); }

    inline bool GetIntegrality(Index iv) const { return domain_.GetIntegrality(iv); }
    inline void SetIntegrality(Index iv, bool b) { domain_.SetIntegrality(iv, b); }

    bool IsIntegerFeasible(const std::vector<Scalar>& x);
    bool IsFeasible(const std::vector<Scalar>& x);

    void AddSlacks();
    void PruneSlacks();

    void Resize(Index m, Index n) {
        obj_.coefficients.resize(n);
        Ac_.Resize(m, n);
        Ar_.Resize(m, n);
        rhs_.resize(m);
        domain_.Resize(n);
    }

    void PrepareConstraint(const SparseVector& lhs, const Bounds& rhs) {
        // Unsafe, need to finalize after this;
        this->Ar_.GetRows().push_back(lhs);
        this->rhs_.push_back(rhs);
    }
    void FinalizeAc();

    Solution PrepareSolution(const LpStatus status, const std::vector<Scalar>& x) const;

    friend std::ostream& operator<<(std::ostream& os, const MilpModel& model);

   private:
    Objective obj_;
    SparseColMatrix Ac_;
    SparseRowMatrix Ar_;
    std::vector<Bounds> rhs_;
    Domain domain_;

    bool has_slacks_ = false;
};

}  // namespace reshala
