#pragma once

#include "Matrix/matrix.h"

template<class T>
Matrix<T> GetRotationMatrix(const Matrix<T>& a) {
  if (!a.IsVector() || std::max(a.Rows(), a.Cols()) != 2) {
    throw std::runtime_error("Matrix a with size " + PairToString(a.Size())
                                 + " should be 1x2 or 2x1 vector");
  }
  Matrix<T> ans(2);
  auto sqrt = std::sqrt(a(0) * a(0) + a(1) * a(1));
  auto cos = a(0) / sqrt;
  auto sin = a(1) / sqrt;
  ans(0, 0) = ans(1, 1) = cos;
  ans(1, 0) = -sin;
  ans(0, 1) = sin;
  return ans;
}
