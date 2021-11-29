#pragma once

#include "Matrix/matrix.h"

template<class T>
Matrix<T> FrobeniusForm(Matrix<T> a) {
  if (!a.IsSquare()) {
    throw std::invalid_argument(
        "Matrix of size " + PairToString(a.Size()) + " is not square.");
  }
  auto n = a.Rows();
  std::vector<int> reindex(n);
  std::iota(reindex.begin(), reindex.end(), 0);

  for (int i = n - 1; i > 0; i--) {
    int row = reindex[i];
    int col = reindex[i - 1];
    int index_of_max = 0;
    if (i == n - 1) {
      index_of_max = n - 1;
    }
    for (int j = 0; j < i; j++) {
      auto index = reindex[j];
      if (std::abs(a.Row(row)(reindex[index_of_max]))
          < std::abs(a.Row(row)(index))) {
        index_of_max = j;
      }
    }
    std::swap(reindex[i - 1], reindex[index_of_max]);
    row = reindex[i];
    col = reindex[i - 1];
    if (std::abs(a(row, col)) < Matrix<T>::GetEps()) {
      continue;
    }
    auto d = a.Col(col)(row);
    a.Col(col).SubMatrix(0, 0, i + 1, -1) /= d;
    a.Row(col) *= d;
    for (int j = 0; j < n; j++) {
      auto index = reindex[j];
      if (col == index) {
        continue;
      }
      auto d = a.Col(index)(row);
      a.Col(index).SubMatrix(0, 0, i + 1, -1) -=
          a.Col(col).SubMatrix(0, 0, i + 1, -1) * d;
      a.Row(col) += a.Row(index) * d;
    }
    Matrix<T> aa(a.Rows(), a.Cols());
    for (int k = 0; k < n; k++) {
      for (int j = 0; j < n; j++) {
        aa(k, j) = a(reindex[k], reindex[j]);
      }
    }
    std::cerr << aa;
    // break;
  }
  Matrix<T> aa(a.Rows(), a.Cols());
  for (int k = 0; k < n; k++) {
    for (int j = 0; j < n; j++) {
      aa(k, j) = a(reindex[k], reindex[j]);
    }
  }
  return aa;
}
