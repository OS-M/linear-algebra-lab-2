#pragma once

#include "Matrix/matrix.h"

template<class T>
std::pair<Matrix<T>, Matrix<T>> QrDecompose(Matrix<T> a, Matrix<T> b) {
  if (!b.IsColVector()) {
    throw std::invalid_argument(
        "Matrix of size " + PairToString(b.Size()) + " is not col.");
  }
  if (b.Rows() != a.Rows()) {
    throw std::invalid_argument(
        "Bad matrix sizes: " + PairToString(a.Size()) + " and "
            + PairToString(b.Size()));
  }
  for (int i = 0; i < a.Cols(); i++) {
    auto col_vector = a.SubMatrix(i, i, -1, 1);
    auto old_col_vector = col_vector;
    col_vector(0) = (col_vector(0) < 0 ? 1 : -1) * EuclideanNorm<T>(col_vector);
    col_vector.SubMatrix(1, 0, -1, 1) *= 0;
    auto norm = EuclideanNorm<T>(old_col_vector - col_vector);
    if (norm < Matrix<T>::GetEps()) {
      continue;
    }
    auto w = (old_col_vector - col_vector) / norm;

    for (int j = i + 1; j < a.Cols(); j++) {  // for columns
      auto cur_col_vector = a.SubMatrix(i, j, -1, 1);
      cur_col_vector -= 2 * cur_col_vector.ScalarProduct(w) * w;
    }
    b.SubMatrix(i, 0, -1, -1) -=
        2 * b.SubMatrix(i, 0, -1, -1).ScalarProduct(w) * w;
  }
  return {a, b};
}
