#pragma once

#include "reshala/lina/basis.h"
#include "reshala/lina/core/dense_matrix.h"
#include "reshala/lina/core/sparse_matrix.h"

namespace reshala {

struct Eta {
    SparseVector eta;
    Index ir;
    explicit Eta(SparseVector sv, Index i) : eta{std::move(sv)}, ir{i} {}
};

class Lina {
   public:
    Lina() {}
    Lina(const SparseColMatrix* Ac, const SparseRowMatrix* Ar, const LpBasis* basis)
        : Ac_(Ac),
          Ar_(Ar),
          m(Ac->GetNRows()),
          n(Ac->GetNCols()),
          basis_(basis),
          Binv_(m, m),
          Lr(m, m),
          Ur(m, m),
          Lc(m, m),
          Uc(m, m) {
        Init();
    }
    Lina& operator=(const Lina&) = default;
    void Init();

    void Btran(Index iv, DenseVector& res);
    void Ftran(Index iv, DenseVector& res);
    void Update(Index iv_leaving, Index iv_entering);

   private:
    enum class UpdType { kDluSm, kSluSm, kSlu };
    static const UpdType ut = UpdType::kSlu;
    static const Index kMaxUpdates = 10;

    const SparseColMatrix* Ac_;
    const SparseRowMatrix* Ar_;
    Index m;
    Index n;

    const LpBasis* basis_;
    DenseMatrix Binv_;

    Index n_updates_;

    // Dense
    bool InvertD();
    bool BuildBinv();  // From sparse LU
    void SherMor(Index iv_leaving, Index iv_entering);
    void BtranD(Index iv, DenseVector& res);
    void FtranD(Index iv, DenseVector& res);

    // Sparse
    SparseRowMatrix Lr, Ur;
    SparseColMatrix Lc, Uc;
    std::vector<Index> row_perm;
    std::vector<Index> row_perm_inv;
    std::vector<Eta> etas;

    bool SparseLU(const SparseRowMatrix& A);
    bool InvertS();

    void BtranS(Index iv, DenseVector& res);
    void SolveUt(DenseVector& x);
    void SolveLt(DenseVector& x);
    void FtranS(Index iv, DenseVector& res);
    void SolveL(DenseVector& x);
    void SolveU(DenseVector& x);
};

}  // namespace reshala
