#pragma once

#include "Matrix/matrix.h"
#include "gauss.h"
#include "lu_decompose.h"

template<class T>
std::pair<std::complex<T>, std::complex<T>> SolveQuadraticEquation(
    T a, T b, T c) {
  auto discrim = std::complex<T>(b * b - 4 * a * c);
  auto discrim_sqrt = std::sqrt(discrim);
  return {(-b - discrim_sqrt) / (2. * a), (-b + discrim_sqrt) / (2. * a)};
}

template<class T>
std::pair<std::complex<T>,
          std::complex<T>> ExtractEigenvalues2x2(const Matrix<T>& a) {
  if (a.Size() != std::make_pair(2, 2)) {
    throw std::invalid_argument(
        "Matrix a of size " + PairToString(a.Size()) + " should be (2;2)");
  }
  auto b = -a(0, 0) - a(1, 1);
  auto c = a(0, 0) * a(1, 1) - a(0, 1) * a(1, 0);
  return SolveQuadraticEquation(b / b, b, c);
}

template<class T>
std::vector<Matrix<std::complex<T>>> FindEigenvectorsByValues(
    const Matrix<std::complex<T>>& matrix,
    const std::vector<std::complex<T>>& values) {
  Matrix<std::complex<T>> a(matrix.Rows(), matrix.Cols());
  Matrix<std::complex<T>> b(matrix.Rows(), 1);
  std::vector<Matrix<std::complex<T>>> ans;

  // auto[l, u] = LuDecompose(a);

  for (auto value: values) {
    for (int i = 0; i < matrix.Rows(); i++) {
      for (int j = 0; j < matrix.Cols(); j++) {
        a(i, j) = matrix(i, j);
        if (i == j) {
          a(i, j) -= value;
        }
      }
    }
    auto x = GaussSolve(a, b).first;
    if (std::abs(EuclideanNorm<std::complex<T>>(x)) >
          std::abs(Matrix<std::complex<T>>::GetEps())) {
      ans.push_back(x);
    } else {
      ans.push_back(Matrix<std::complex<T>>(0, 0));
    }
  }
  return ans;
}
