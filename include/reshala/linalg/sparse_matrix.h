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

    Index GetNRows() const { return m_; }
    Index GetNCols() const { return n_; }

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

    Index GetNRows() const { return m_; }
    Index GetNCols() const { return n_; }

    std::vector<SparseVector>& getCols() { return cols_; }
    const std::vector<SparseVector>& getCols() const { return cols_; }

   private:
    Index m_;
    Index n_;
    std::vector<SparseVector> cols_;
};

}  // namespace reshala
