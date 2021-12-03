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
  auto[r1, r2] =
    SolveQuadraticEquation<T>(1, c(1), c(0));

  auto u1 = a * u;
  auto u2 = squared_a * u;

  Matrix<std::complex<T>> v1(u.Rows(), 1);
  Matrix<std::complex<T>> v2(u.Rows(), 1);

  for (int i = 0; i < u.Rows(); i++) {
    v1(i) = u2(i) - r2 * u1(i);
    v2(i) = u1(i) - u2(i) / r1;
  }

  return {{r1, v1}, {r2, v2}};
}

template<class T>
T PowerIterationMethod1Iteration(
    const Matrix<T>& a,
    Matrix<T>& u,
    Matrix<T>& y) {
  y.Assign(a * u);
  u.Assign(y / EuclideanNorm<T>(y));
  auto lambda = (u.ScalarProduct(a * u));
  return lambda;
}

template<class T>
std::pair<bool, Matrix<T>> PowerIterationMethod1IterationConverges(
    const Matrix<T>& a,
    int iters = 10) {
  int n = a.Rows();
  auto y = Matrix<T>(n, 1);
  y(0) = 1;
  auto u = y / EuclideanNorm<T>(y);
  auto lambda = u.ScalarProduct(a * u);
  auto prev_lambda = lambda;
  auto prev_diff = 1e18;
  for (int i = 0; i < iters; i++) {
    lambda = __internal::PowerIterationMethod1Iteration(a, u, y);
    if (std::abs(prev_lambda - lambda) <= prev_diff) {
      prev_diff = std::abs(prev_lambda - lambda);
      prev_lambda = lambda;
      continue;
    } else {
      return {false, y};
    }
  }
  return {true, y};
}

template<class T>
T PowerIterationMethod2Iteration(
    const Matrix<T>& a,
    const Matrix<T>& a2,
    Matrix<T>& u,
    Matrix<T>& y) {
  y.Assign(a2 * u);
  u.Assign(y / EuclideanNorm<T>(y));
  auto lambda = (u.ScalarProduct(a2 * u));
  return lambda;
}

template<class T>
std::pair<bool, Matrix<T>> PowerIterationMethod2IterationConverges(
    const Matrix<T>& a,
    const Matrix<T>& a2,
    int iters = 10) {
  int n = a.Rows();
  auto y = Matrix<T>(n, 1);
  y(0) = 1;
  auto u = y / EuclideanNorm<T>(y);
  auto lambda = u.ScalarProduct(a * u);
  auto prev_lambda = lambda;
  auto prev_diff = 1e18;
  for (int i = 0; i < iters; i++) {
    lambda = __internal::PowerIterationMethod2Iteration(a, a2, u, y);
    if (lambda > 0 && std::abs(prev_lambda - lambda) <= prev_diff) {
      prev_diff = std::abs(prev_lambda - lambda);
      prev_lambda = lambda;
      continue;
    } else {
      return {false, y};
    }
  }
  return {true, y};
}

template<class T>
std::pair<T, Matrix<T>> PowerMethodEigenvalues1(
    const Matrix<T>& a,
    Matrix<T> y,
    int* iters = nullptr,
    int max_iters = 100) {
  if (!a.IsSquare()) {
    throw std::invalid_argument(
        "Matrix of size " + PairToString(a.Size()) + " is not square.");
  }
  int n = a.Rows();
  y(0) = 1;
  auto u = y / EuclideanNorm<T>(y);
  auto lambda = u.ScalarProduct(a * u);
  int iter = 0;
  while (EuclideanNorm<T>(a * u - lambda * u) > Matrix<T>::GetEps()) {
    lambda = __internal::PowerIterationMethod1Iteration(a, u, y);
    iter++;
    if (iter > max_iters) {
      break;
    }
  }
  if (iters) {
    *iters = iter + 1;
    if (iter >= max_iters) {
      *iters = -1;
    }
  }
  return {lambda, u};
}

template<class T>
std::vector<std::pair<std::complex<T>,
                      Matrix<std::complex<T>>>> PowerMethodEigenvalues2(
    const Matrix<T>& a,
    Matrix<T> y,
    int* iters = nullptr,
    int max_iters = 100) {
  if (!a.IsSquare()) {
    throw std::invalid_argument(
        "Matrix of size " + PairToString(a.Size()) + " is not square.");
  }
  auto a2 = a * a;
  int n = a.Rows();
  auto u = y / EuclideanNorm<T>(y);
  auto lambda = u.ScalarProduct(a * u);
  auto v1 = u;
  auto v2 = u;
  int iter = 0;
  while (
      EuclideanNorm<T>(a * v1 - std::sqrt(lambda) * v1) > Matrix<T>::GetEps() ||
      EuclideanNorm<T>(a * v2 + std::sqrt(lambda) * v2) > Matrix<T>::GetEps()) {
    lambda = __internal::PowerIterationMethod2Iteration(a, a2, u, y);
    v1 = (std::sqrt(lambda) * a * u + a2 * u) / (2 * lambda);
    v2 = (-std::sqrt(lambda) * a * u + a2 * u) / (2 * lambda);
    iter++;
    if (iter > max_iters) {
      break;
    }
  }
  if (iters) {
    *iters = iter + 1;
    if (iter >= max_iters) {
      *iters = -1;
    }
  }

  return {std::make_pair(std::sqrt(lambda), v1.ToComplex()),
          std::make_pair(-std::sqrt(lambda), v2.ToComplex())};
}

template<class T>
std::vector<std::pair<std::complex<T>,
                      Matrix<std::complex<T>>>> PowerMethodEigenvalues3(
    const Matrix<T>& a,
    int* iters = nullptr,
    int max_iters = 100) {
  if (!a.IsSquare()) {
    throw std::invalid_argument(
        "Matrix of size " + PairToString(a.Size()) + " is not square.");
  }
  auto squared_a = a * a;
  int n = a.Rows();

  auto complex_a = a.ToComplex();
  auto complex_squared_a = squared_a.ToComplex();

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

  while (error1 > std::norm(Matrix<T>::GetEps()) ||
      error2 > std::norm(Matrix<T>::GetEps())) {
    auto[p1, p2] =
    __internal::PowerMethodEigenvaluesComplexIteration(
        complex_a, complex_squared_a, complex_y, complex_u);
    r1 = p1.first;
    v1 = p1.second;
    r2 = p2.first;
    v2 = p2.second;
    error1 = std::norm(
        EuclideanNorm<std::complex<T>>(complex_a * v1 - r1 * v1));
    error2 = std::norm(
        EuclideanNorm<std::complex<T>>(complex_a * v2 - r2 * v2));
    iter++;
    if (iter > max_iters) {
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

  // std::cerr << error1 << ' ' << error2 << '\n';

  std::vector<std::pair<std::complex<T>,
                        Matrix<std::complex<T>>>> ans;
  if (error1 < std::norm(Matrix<T>::GetEps())
      && std::norm(EuclideanNorm<std::complex<T>>(v1))
          > std::norm(Matrix<std::complex<T>>::GetEps())) {
    ans.emplace_back(r1, v1);
  }
  if (error2 < std::norm(Matrix<T>::GetEps())
      && std::norm(EuclideanNorm<std::complex<T>>(v2))
          > std::norm(Matrix<std::complex<T>>::GetEps())) {
    ans.emplace_back(r2, v2);
  }
  return ans;
}

}

template<class T>
std::vector<std::pair<std::complex<T>,
                      Matrix<std::complex<T>>>> PowerMethodEigenvalues(
    const Matrix<T>& a,
    int* iters = nullptr,
    int max_iters = 100) {
  {
    auto[is_converging, y] =
      __internal::PowerIterationMethod1IterationConverges(a);
    if (is_converging) {
      auto[e, v] =
        __internal::PowerMethodEigenvalues1(a, y, iters, max_iters);
      std::cerr << "Method 1\n";
      return {std::make_pair(e, v.ToComplex())};
    }
  }
  {
    auto[is_converging, y] =
      __internal::PowerIterationMethod2IterationConverges(a, a * a);
    if (is_converging) {
      std::cerr << "Method 2\n";
      return __internal::PowerMethodEigenvalues2(a, y, iters, max_iters);
    }
  }
  std::cerr << "Method 3\n";
  return __internal::PowerMethodEigenvalues3(a, iters, max_iters);
}
