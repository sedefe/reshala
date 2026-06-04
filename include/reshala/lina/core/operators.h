#pragma once

#include "dense_matrix.h"
#include "sparse_matrix.h"

namespace reshala {

void dot(const Scalar* ar1, const SparseVector& sv2, Scalar& res);
void dot(const SparseVector& sv1, const SparseVector& sv2, Scalar& res);
void dot(const DenseVector& dv1, const SparseVector& sv2, Scalar& res);
void dot(const DenseVector& dv1, const DenseVector& dv2, Scalar& res);

void MulScmDv(const SparseColMatrix& scm, const DenseVector& sv, DenseVector& res);
void MulDmSv(const DenseMatrix& dm, const SparseVector& sv, DenseVector& res);

bool Invert(DenseMatrix& dm);

}  // namespace reshala
