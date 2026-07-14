#pragma once

#include "reshala/lina/basis.h"
#include "reshala/lina/core/sparse_matrix.h"

namespace reshala {

struct Eta {
    SparseVector eta;
    Index i_col;
    Scalar diag;
    Eta(SparseVector sv, Index i) : eta{std::move(sv)}, i_col{i} { diag = eta.At(i_col); }
};

struct LinaStats {
    Index n_lus = 0;
    Index n_updates = 0;

    Index total_nnz_b = 0;
    Index total_nnz_l = 0;
    Index total_nnz_u = 0;
};
std::ostream& operator<<(std::ostream& os, const LinaStats& stats);

class Lina {
   public:
    Lina() : ftran_res(0) {}
    Lina(const SparseColMatrix& Ac, const SparseRowMatrix& Ar, const LpBasis* basis);
    Lina& operator=(const Lina&) = default;

    bool Refactor();
    inline Index GetAge() const { return etas.size(); }

    void Btran(const SparseVector& b, DenseVector& res) const;
    void Btran(DenseVector& x, DenseVector& res) const;   // NB: x is modified!
    void Ftran(const SparseVector& b, DenseVector& res);  // NB: stores the result
    void Ftran(const DenseVector& x, DenseVector& res) const;
    void Update(Index iv_leaving, Index iv_entering);

    inline const LinaStats& GetStats() const { return stats; }

   private:
    enum class UpdType { kSlu, kSluPf };
    static const UpdType ut = UpdType::kSluPf;

    LinaStats stats;

    SparseColMatrix Ac_;
    SparseRowMatrix Ar_;
    Index m;
    Index n;

    const LpBasis* basis_;

    SparseRowMatrix Lr, Ur;
    SparseColMatrix Lc, Uc;
    DenseVector u_diag;

    std::vector<Index> row_perm;
    std::vector<Index> row_perm_inv;
    std::vector<Eta> etas;

    SparseVector ftran_res;

    void ProdForm(Index iv_leaving, Index iv_entering);

    void SolveUt(DenseVector& x) const;
    void SolveLt(DenseVector& x) const;
    void EtaBtran(const Eta& eta, DenseVector& x) const;

    void SolveL(DenseVector& x) const;
    void SolveU(DenseVector& x) const;
    void EtaFtran(const Eta& eta, DenseVector& x) const;
};

}  // namespace reshala
