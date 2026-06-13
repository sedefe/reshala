#pragma once

#include <algorithm>
#include <ostream>
#include <vector>

#include "dense_vector.h"
#include "sparse_vector.h"

namespace reshala {

class SparseRowMatrix {
   public:
    SparseRowMatrix(Index m = 0, Index n = 0) : m_(m), n_(n), rows_(m, SparseVector(n)) {}

    Index GetNRows() const { return m_; }
    Index GetNCols() const { return n_; }
    Index GetNnz() const {
        Index nnz = 0;
        for (const auto& row : GetRows()) {
            nnz += row.Size();
        }
        return nnz;
    }

    void Resize(Index m, Index n) {
        if (n != n_) {
            for (auto& row : rows_) row.SetDim(n);
        }
        m_ = m;
        n_ = n;
        rows_.resize(m, n);
    }
    inline void Clear() {
        for (auto& row : rows_) row.Clear();
    }

    std::vector<SparseVector>& GetRows() { return rows_; }
    const std::vector<SparseVector>& GetRows() const { return rows_; }

    SparseVector& GetRow(Index i) { return rows_[i]; }
    const SparseVector& GetRow(Index i) const { return rows_[i]; }

   private:
    Index m_;
    Index n_;
    std::vector<SparseVector> rows_;
};

class SparseColMatrix {
   public:
    SparseColMatrix(Index m = 0, Index n = 0) : m_(m), n_(n), cols_(n, SparseVector(m)) {}

    Index GetNRows() const { return m_; }
    Index GetNCols() const { return n_; }
    Index GetNnz() const {
        Index nnz = 0;
        for (const auto& col : GetCols()) {
            nnz += col.Size();
        }
        return nnz;
    }

    void Resize(Index m, Index n) {
        if (m != m_) {
            for (auto& col : cols_) col.SetDim(m);
        }
        m_ = m;
        n_ = n;
        cols_.resize(n, m);
    }
    inline void Clear() {
        for (auto& col : cols_) col.Clear();
    }

    std::vector<SparseVector>& GetCols() { return cols_; }
    const std::vector<SparseVector>& GetCols() const { return cols_; }

    SparseVector& GetCol(Index i) { return cols_[i]; }
    const SparseVector& GetCol(Index i) const { return cols_[i]; }

   private:
    Index m_;
    Index n_;
    std::vector<SparseVector> cols_;
};

std::ostream& operator<<(std::ostream& os, const SparseRowMatrix& svm);

void Srm2Scm(const SparseRowMatrix& srm, SparseColMatrix& scm);

}  // namespace reshala
