#pragma once

#include "Matrix/matrix.h"
#include "euclidean_norm.h"

template<class T>
Matrix<T> ReflectionsHessenberg(Matrix<T> a) {
  if (!a.IsSquare()) {
    throw std::invalid_argument(
        "Matrix of size " + PairToString(a.Size()) + " is not square.");
  }
  int n = a.Rows();
  for (int i = 0; i < n - 2; i++) {
    auto col_vector = a.SubMatrix(i + 1, i, -1, 1);
    auto old_col_vector = col_vector;
    col_vector(0) = (col_vector(0) < 0 ? 1 : -1) * EuclideanNorm<T>(col_vector);
    col_vector.SubMatrix(1, 0, -1, 1) *= 0;
    auto norm = EuclideanNorm<T>(old_col_vector - col_vector);
    if (norm < Matrix<T>::GetEps()) {
      continue;
    }
    auto w = (old_col_vector - col_vector) / norm;

    for (int j = i + 1; j < n; j++) {  // for columns
      auto cur_col_vector = a.SubMatrix(i + 1, j, -1, 1);
      cur_col_vector -= 2 * cur_col_vector.ScalarProduct(w) * w;
    }
    w = w.Transposed();
    for (int j = 0; j < n; j++) {  // for rows
      auto cur_row_vector = a.SubMatrix(j, i + 1, 1, -1);
      cur_row_vector -= 2 * cur_row_vector.ScalarProduct(w) * w;
    }
  }
  return a;
}
