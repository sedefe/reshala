#pragma once

#include "dense_matrix.h"
#include "sparse_matrix.h"

namespace reshala {

SparseVector sub(const SparseVector& x, const SparseVector& y);
void dot(const Scalar* ar1, const SparseVector& sv2, Scalar& res);
void dot(const SparseVector& sv1, const SparseVector& sv2, Scalar& res);
void dot(const DenseVector& dv1, const SparseVector& sv2, Scalar& res);
void dot(const DenseVector& dv1, const DenseVector& dv2, Scalar& res);

void MulScmSv(const SparseColMatrix& scm, const SparseVector& sv, DenseVector& res);
void MulScmDv(const SparseColMatrix& scm, const DenseVector& sv, DenseVector& res);
void MulDvDm(const DenseVector& dv, const DenseMatrix& dm, DenseVector& res);
void MulDmDv(const DenseMatrix& dm, const DenseVector& dv, DenseVector& res);
void MulDmSv(const DenseMatrix& dm, const SparseVector& sv, DenseVector& res);

}  // namespace reshala
