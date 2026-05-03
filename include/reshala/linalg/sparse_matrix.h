#pragma once

#include <algorithm>
#include <cassert>
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
        for (const auto& row : getRows()) {
            nnz += row.size();
        }
        return nnz;
    }

    void resize(Index m, Index n) {
        m_ = m;
        n_ = n;
        rows_.resize(m, n);
    }

    std::vector<SparseVector>& getRows() { return rows_; }
    const std::vector<SparseVector>& getRows() const { return rows_; }

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
        for (const auto& col : getCols()) {
            nnz += col.size();
        }
        return nnz;
    }

    void resize(Index m, Index n) {
        m_ = m;
        n_ = n;
        cols_.resize(n, m);
    }

    std::vector<SparseVector>& getCols() { return cols_; }
    const std::vector<SparseVector>& getCols() const { return cols_; }

   private:
    Index m_;
    Index n_;
    std::vector<SparseVector> cols_;
};

void Srm2Scm(const SparseRowMatrix& srm, SparseColMatrix& scm);

}  // namespace reshala
