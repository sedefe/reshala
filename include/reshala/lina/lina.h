#pragma once

#include "reshala/lina/basis.h"
#include "reshala/lina/core/dense_matrix.h"
#include "reshala/lina/core/sparse_matrix.h"
#include "reshala/lina/scaling.h"

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
};
std::ostream& operator<<(std::ostream& os, const LinaStats& stats);

class Lina {
   public:
    Lina() : ftran_res(0) {}
    Lina(const SparseColMatrix& Ac, const SparseRowMatrix& Ar, const LpBasis* basis);
    Lina& operator=(const Lina&) = default;
    void Scale();

    bool Refactor();
    inline Index GetAge() const { return etas.size(); }

    void Btran(Index iv, DenseVector& res);
    void Ftran(Index iv, DenseVector& res);
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
    Scaling scaling;

    const LpBasis* basis_;

    SparseRowMatrix Lr, Ur;
    SparseColMatrix Lc, Uc;
    DenseVector u_diag;

    std::vector<Index> row_perm;
    std::vector<Index> row_perm_inv;
    std::vector<Eta> etas;

    SparseVector ftran_res;

    void ProdForm(Index iv_leaving, Index iv_entering);

    void SolveUt(DenseVector& x);
    void SolveLt(DenseVector& x);
    void EtaBtran(const Eta& eta, DenseVector& x);

    void SolveL(DenseVector& x);
    void SolveU(DenseVector& x);
    void EtaFtran(const Eta& eta, DenseVector& x);
};

}  // namespace reshala
