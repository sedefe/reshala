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

    Index getNRows() const { return m_; }
    Index getNCols() const { return n_; }
    Index getNnz() {
        Index nnz = 0;
        for (const auto& row : GetRows()) {
            nnz += row.Size();
        }
        return nnz;
    }

    void Resize(Index m, Index n) {
        m_ = m;
        n_ = n;
        rows_.resize(m, n);
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

    Index getNRows() const { return m_; }
    Index getNCols() const { return n_; }
    Index getNnz() {
        Index nnz = 0;
        for (const auto& col : GetCols()) {
            nnz += col.Size();
        }
        return nnz;
    }

    void Resize(Index m, Index n) {
        m_ = m;
        n_ = n;
        cols_.resize(n, m);
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

void Srm2Scm(const SparseRowMatrix& srm, SparseColMatrix& scm);

}  // namespace reshala
