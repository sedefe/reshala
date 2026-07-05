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
    void Reset();
    void Scale();

    void Btran(Index iv, DenseVector& res);
    void Ftran(Index iv, DenseVector& res);
    void Update(Index iv_leaving, Index iv_entering);

    inline const LinaStats& GetStats() const { return stats; }

   private:
    enum class UpdType { kDluSm, kSluSm, kSlu, kSluPf };
    static const UpdType ut = UpdType::kSluPf;
    static const Index kMaxUpdates = 50;

    LinaStats stats;

    SparseColMatrix Ac_;
    SparseRowMatrix Ar_;
    Index m;
    Index n;
    Scaling scaling;

    const LpBasis* basis_;
    DenseMatrix Binv_;

    Index n_updates_ = 0;

    // Dense
    bool InvertD();
    bool SparseLu2Binv();  // From sparse LU
    void SherMor(Index iv_leaving, Index iv_entering);
    void BtranD(Index iv, DenseVector& res);
    void FtranD(Index iv, DenseVector& res);

    // Sparse
    SparseRowMatrix Lr, Ur;
    SparseColMatrix Lc, Uc;

    std::vector<Index> row_perm;
    std::vector<Index> row_perm_inv;
    std::vector<Eta> etas;

    SparseVector ftran_res;

    bool SparseLU();
    void ProdForm(Index iv_leaving, Index iv_entering);

    void BtranS(Index iv, DenseVector& res);
    void SolveUt(DenseVector& x);
    void SolveLt(DenseVector& x);
    void EtaBtran(const Eta& eta, DenseVector& x);

    void FtranS(Index iv, DenseVector& res);
    void SolveL(DenseVector& x);
    void SolveU(DenseVector& x);
    void EtaFtran(const Eta& eta, DenseVector& x);
};

}  // namespace reshala
