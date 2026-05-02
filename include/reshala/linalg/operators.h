#pragma once

#include "sparse_matrix.h"

namespace reshala {

void dot(const SparseVector& sv1, const SparseVector& sv2, Scalar& res);
void dot(const DenseVector& dv1, const SparseVector& sv2, Scalar& res);
void dot(const DenseVector& dv1, const DenseVector& dv2, Scalar& res);

void MulScmSv(const SparseColMatrix& scm, const SparseVector& sv, DenseVector& res);

}  // namespace reshala
