#pragma once

#include "reshala/linalg/sparse_matrix.h"
#include "reshala/model/domain.h"
#include "reshala/model/objective.h"
#include "reshala/model/solution.h"
#include "reshala/model/utils.h"

namespace reshala {

class MilpModel {
   public:
    MilpModel() {}

    Index GetNCons() const { return Ac_.GetNRows(); }
    Index GetNVars() const { return Ac_.GetNCols(); }
    Index GetNnz() const;

    const Objective& GetObj() const { return obj_; }
    Objective& GetObj() { return obj_; }

    inline const SparseColMatrix& GetAc() const { return Ac_; }
    inline SparseColMatrix& GetAc() { return Ac_; }
    inline const SparseVector& GetCol(Index i) const { return Ac_.GetCol(i); }
    inline SparseVector& GetCol(Index i) { return Ac_.GetCol(i); }

    inline const SparseRowMatrix& GetAr() const { return Ar_; }
    inline SparseRowMatrix& GetAr() { return Ar_; }
    inline const SparseVector& GetRow(Index i) const { return Ar_.GetRow(i); }
    inline SparseVector& GetRow(Index i) { return Ar_.GetRow(i); }

    const std::vector<Bounds>& GetRhs() const { return rhs_; }
    std::vector<Bounds>& GetRhs() { return rhs_; }
    const Bounds& GetRhs(Index ic) const { return rhs_[ic]; }
    Bounds& GetRhs(Index ic) { return rhs_[ic]; }

    const Domain& GetVars() const { return domain_; }
    Domain& GetVars() { return domain_; }

    inline const Domain& GetDomain() const { return domain_; }
    inline Domain& GetDomain() { return domain_; }
    inline void SetDomain(const Domain& domain) { domain_ = domain; }
    inline const Bounds& GetBounds(Index iv) const { return domain_.GetBounds(iv); }
    inline void SetBounds(Index iv, const Bounds& bnd) { domain_.SetBounds(iv, bnd); }
    inline const BndType& GetType(Index iv) const { return domain_.GetType(iv); }

    inline bool GetIntegrality(Index iv) const { return domain_.GetIntegrality(iv); }
    inline void SetIntegrality(Index iv, bool b) { domain_.SetIntegrality(iv, b); }
    inline bool IsBinary(Index iv) const {
        return GetIntegrality(iv) && (domain_.GetBounds(iv).ri == 1) &&
               (domain_.GetBounds(iv).le == 0);
    }

    bool IsIntegerFeasible(const std::vector<Scalar>& x) const;
    FeasibilityReport GetFeasReport(const std::vector<Scalar>& x) const;

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
