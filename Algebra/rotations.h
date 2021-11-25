#pragma once

#include "Matrix/matrix.h"

template<class T>
std::pair<T, T> GetRotationMatrix(const Matrix<T>& a) {
  if (!a.IsVector() || std::max(a.Rows(), a.Cols()) != 2) {
    throw std::runtime_error("Matrix a with size " + PairToString(a.Size())
                                 + " should be 1x2 or 2x1 vector");
  }
  Matrix<T> ans(2);
  auto sqrt = std::sqrt(a(0) * a(0) + a(1) * a(1));
  auto cos = a(0) / sqrt;
  auto sin = a(1) / sqrt;
  return {sin, cos};
}

template<class T>
void ApplyRotation(Matrix<T>& a, T sin, T cos, int iter) {
  if (!a.IsSquare()) {
    throw std::invalid_argument(
        "Matrix of size " + PairToString(a.Size()) + " is not square.");
  }
  auto n = a.Rows();
  for (int i = iter; i < n; i++) {
    T e1 = a(iter, i) * cos + a(iter + 1, i) * sin;
    T e2 = a(iter, i) * (-sin) + a(iter + 1, i) * cos;
    a(iter, i) = e1;
    a(iter + 1, i) = e2;
  }
}

template<class T>
void ApplyTransposedRotation(Matrix<T>& a, T sin, T cos, int iter) {
  if (!a.IsSquare()) {
    throw std::invalid_argument(
        "Matrix of size " + PairToString(a.Size()) + " is not square.");
  }
  for (int i = 0; i < iter + 2; i++) {
    T e1 = a(i, iter) * cos + a(i, iter + 1) * sin;
    T e2 = a(i, iter) * (-sin) + a(i, iter + 1) * cos;
    a(i, iter) = e1;
    a(i, iter + 1) = e2;
  }
}
