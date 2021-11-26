#pragma once

#include "Matrix/matrix.h"
#include "euclidean_norm.h"

template<class T>
std::vector<std::pair<T, Matrix<T>>> PowerMethodEigenvalues(
    const Matrix<T>& a) {
  if (!a.IsSquare()) {
    throw std::invalid_argument(
        "Matrix of size " + PairToString(a.Size()) + " is not square.");
  }
  auto a2 = a * a;
  int n = a.Rows();
  auto y = Matrix<T>(n, 1);
  y(0) = 1;
  auto u = y / EuclideanNorm<T>(y);
  auto lambda = u.ScalarProduct(a * u);
  // while (EuclideanNorm<T>(a * u - lambda * u) > Matrix<T>::eps) {
  for (int i = 0; i < 10; i++) {
    y.Assign(a2 * u);
    u.Assign(y / EuclideanNorm<T>(y));
    lambda = (u.ScalarProduct(a2 * u));
  }
  std::cerr << std::sqrt(lambda) << '\n' << u;
  auto aa = (std::sqrt(lambda) * a * u + a2 * u) / (2 * lambda);
  auto aa2 = (-std::sqrt(lambda) * a * u + a2 * u) / (2 * lambda);
  std::cerr << aa / EuclideanNorm<T>(aa);
  std::cerr << aa2 / EuclideanNorm<T>(aa2);
  return {};
}
