#pragma once

#include "Matrix/matrix.h"
#include "euclidean_norm.h"
#include "minimal_square_problem.h"
#include "eigenvalues.h"

namespace __internal {

template<class T>
std::pair<std::pair<std::complex<T>, Matrix<std::complex<T>>>,
          std::pair<std::complex<T>, Matrix<std::complex<T>>>>
PowerMethodEigenvaluesComplexIteration(
    const Matrix<std::complex<T>>& a,
    const Matrix<std::complex<T>>& squared_a,
    Matrix<std::complex<T>>& y,
    Matrix<std::complex<T>>& u) {
  y.Assign(a * u);
  u.Assign(y / EuclideanNorm<std::complex<T>>(y));

  Matrix<T> l(u.Rows(), 2);
  auto au = a * u;
  for (int i = 0; i < u.Rows(); i++) {
    l(i, 0) = u(i).real();
    l(i, 1) = au(i).real();
  }
  auto rc = -1 * squared_a * u;
  Matrix<T> r(u.Rows(), 1);
  for (int i = 0; i < u.Rows(); i++) {
    r(i) = rc(i).real();
  }
  auto c = MinimalSquareProblem(l, r);
  auto[r1, r2] = SolveQuadraticEquation<T>(1, c(1), c(0));

  auto u1 = a * u;
  auto u2 = squared_a * u;

  Matrix<std::complex<T>> v1(u.Rows(), 1);
  Matrix<std::complex<T>> v2(u.Rows(), 1);

  for (int i = 0; i < u.Rows(); i++) {
    v1(i) = u2(i) - r2 * u1(i);
    v2(i) = u1(i) - u2(i) / r1;
  }

  // auto v1 = (r1 * a * u + squared_a * u) / (std::complex<T>(2) * r1 * r1);
  // auto v2 = (r2 * a * u + squared_a * u) / (std::complex<T>(2) * r2 * r2);
  return {{r1, v1}, {r2, v2}};
}

}

template<class T>
std::vector<std::pair<std::complex<T>,
                      Matrix<std::complex<T>>>> PowerMethodEigenvalues(
    const Matrix<T>& a,
    int* iters = nullptr,
    int max_iters = 100) {
  if (!a.IsSquare()) {
    throw std::invalid_argument(
        "Matrix of size " + PairToString(a.Size()) + " is not square.");
  }
  auto squared_a = a * a;
  int n = a.Rows();

  Matrix<std::complex<T>> complex_a(a.Rows(), a.Cols());
  Matrix<std::complex<T>> complex_squared_a(a.Rows(), a.Cols());
  for (int i = 0; i < a.Rows(); i++) {
    for (int j = 0; j < a.Cols(); j++) {
      complex_a(i, j) = a(i, j);
      complex_squared_a(i, j) = squared_a(i, j);
    }
  }

  Matrix<std::complex<T>> complex_y(n, 1);
  complex_y(0) = 1;
  auto complex_u = complex_y / EuclideanNorm<std::complex<T>>(complex_y);

  T error1 = 1e18;
  T error2 = 1e18;

  std::complex<T> r1;
  std::complex<T> r2;
  Matrix<std::complex<T>> v1;
  Matrix<std::complex<T>> v2;

  int iter = 0;

  while (error1 > Matrix<T>::GetEps() ||
         error2 > Matrix<T>::GetEps()) {
    auto[p1, p2] = __internal::PowerMethodEigenvaluesComplexIteration(
        complex_a, complex_squared_a, complex_y, complex_u);
    r1 = p1.first;
    v1 = p1.second;
    r2 = p2.first;
    v2 = p2.second;
    error1 = std::sqrt(std::norm(
        EuclideanNorm<std::complex<T>>(complex_a * v1 - r1 * v1)));
    error2 = std::sqrt(std::norm(
        EuclideanNorm<std::complex<T>>(complex_a * v2 - r2 * v2)));
    iter++;
    if (iter > max_iters && std::abs(error1 - error2) < Matrix<T>::GetEps()) {
      std::cerr << error1 << ' ' << error2 << '\n';
      break;
    }
  }

  if (iters) {
    *iters = iter + 1;
    if (iter >= max_iters) {
      *iters = -1;
    }
  }

  std::cerr << error1 << ' ' << error2 << '\n';

  std::vector<std::pair<std::complex<T>,
                        Matrix<std::complex<T>>>> ans;
  if (error1 < Matrix<T>::GetEps()) {
    ans.emplace_back(r1, v1);
  }
  if (error2 < Matrix<T>::GetEps()) {
    ans.emplace_back(r2, v2);
  }
  return ans;
}
