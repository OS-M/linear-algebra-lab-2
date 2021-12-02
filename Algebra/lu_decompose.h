#pragma once

#include "Matrix/matrix.h"

template<class T>
std::pair<Matrix<T>, Matrix<T>> LuDecompose(const Matrix<T>& a) {
  if (!a.IsSquare()) {
    throw std::runtime_error("A is not square");
  }

  Matrix<T> l(a.Rows(), a.Cols());
  Matrix<T> u(a.Rows(), a.Cols());

  u = a;
  auto n = a.Rows();
  for (int k = 0; k < n - 1; k++) {
    for (int i = k + 1; i < n; i++) {
      u.At(i, k) /= u.At(k, k);
      for (int j = k + 1; j < n; j++) {
        u.At(i, j) -= u.At(i, k) * u.At(k, j);
      }
    }
  }
  for (int i = 0; i < n; i++) {
    for (int j = 0; j <= i; j++) {
      l.At(i, j) = u.At(i, j);
      if (i == j) {
        l.At(i, j) = 1;
      } else {
        u.At(i, j) = 0;
      }
    }
  }
  return {l, u};
}

template<class T>
Matrix<T> LuSolve(const Matrix<T>& l,
                  const Matrix<T>& u,
                  const Matrix<T>& b) {
  return SolveUxb(u, SolveLxb(l, b));
}
