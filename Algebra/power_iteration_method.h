#pragma once

#include "Matrix/matrix.h"
#include "euclidean_norm.h"
#include "minimal_square_problem.h"
#include "eigenvalues.h"

template<class T>
std::vector<std::pair<T, Matrix<T>>> PowerMethodEigenvalues1(
    const Matrix<T>& a) {
  if (!a.IsSquare()) {
    throw std::invalid_argument(
        "Matrix of size " + PairToString(a.Size()) + " is not square.");
  }
  int n = a.Rows();
  auto y = Matrix<T>(n, 1);
  y(0) = 1;
  auto u = y / EuclideanNorm<T>(y);
  auto lambda = u.ScalarProduct(a * u);
  // while (EuclideanNorm<T>(a * u - lambda * u) > Matrix<T>::GetEps()) {
  for (int i = 0; i < 10; i++) {
    y.Assign(a * u);
    u.Assign(y / EuclideanNorm<T>(y));
    lambda = (u.ScalarProduct(a * u));
  }
  std::cerr << lambda << '\n' << u;
  return {};
}

template<class T>
std::vector<std::pair<T, Matrix<T>>> PowerMethodEigenvalues2(
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
  // while (EuclideanNorm<T>(a * u - lambda * u) > Matrix<T>::GetEps()) {
  for (int i = 0; i < 15; i++) {
    y.Assign(a2 * u);
    u.Assign(y / EuclideanNorm<T>(y));
    lambda = (u.ScalarProduct(a2 * u));
  }
  std::cerr << u;
  auto aa = (std::sqrt(lambda) * a * u + a2 * u) / (2 * lambda);
  auto aa2 = (-std::sqrt(lambda) * a * u + a2 * u) / (2 * lambda);
  std::cerr << std::sqrt(lambda) << '\n';
  std::cerr << aa / EuclideanNorm<T>(aa);
  std::cerr << aa2 / EuclideanNorm<T>(aa2);
  return {};
}

template<class T>
std::vector<std::pair<T, Matrix<T>>> PowerMethodEigenvalues3(
    const Matrix<T>& a) {
  if (!a.IsSquare()) {
    throw std::invalid_argument(
        "Matrix of size " + PairToString(a.Size()) + " is not square.");
  }
  int n = a.Rows();
  auto y = Matrix<T>(n, 1);
  y(0) = 1;
  auto u = y / EuclideanNorm<T>(y);
  auto lambda = u.ScalarProduct(a * u);
  // while (EuclideanNorm<T>(a * u - lambda * u) > Matrix<T>::GetEps()) {
  for (int g = 0; g < 150; g++) {
    y.Assign(a * u);
    u.Assign(y / EuclideanNorm<T>(y));

    Matrix<T> l(u.Rows(), 2);
    l.Col(0).Assign(u);
    l.Col(1).Assign(a * u);
    auto r = -1 * a * a * u;
    auto c = MinimalSquareProblem(l, r);
    auto[r1, r2] = SolveQuadraticEquation<T>(1, c(1), c(0));
    std::cerr << r1 << ' ' << r2 << '\n';
    Matrix<std::complex<T>> a_c(a.Rows(), a.Cols());
    Matrix<std::complex<T>> u_c(u.Rows(), u.Cols());
    for (int i = 0; i < a.Rows(); i++) {
      for (int j = 0; j < a.Cols(); j++) {
        a_c(i, j) = a(i, j);
      }
    }
    for (int i = 0; i < u_c.Rows(); i++) {
      for (int j = 0; j < u_c.Cols(); j++) {
        u_c(i, j) = u(i, j);
      }
    }
    auto aa =
        (r1 * a_c * u_c + a_c * a_c * u_c) / (std::complex<T>(2) * r1 * r1);
    // auto aa2 = (lambda * a * u + a * a * u) / (2 * lambda * lambda);
    std::cerr << norm(EuclideanNorm<std::complex<T>>(a_c * aa - r1 * aa))
              << '\n';
    // std::cerr << aa2 / EuclideanNorm<T>(aa2);
    std::cerr << "--------\n";
  }
  return {};
}

template<class T>
std::vector<Matrix<T>, Matrix<T>> PowerMethodEigenvaluesMethod2Iteration(
    const Matrix<T>& a,
    const Matrix<T>& squared_a,
    Matrix<T>& y,
    Matrix<T>& u,
    T& lambda) {
  y.Assign(squared_a * u);
  u.Assign(y / EuclideanNorm<T>(y));
  lambda = (u.ScalarProduct(squared_a * u));
  return {
      (std::sqrt(lambda) * a * u + squared_a * u) / (2 * lambda),
      (-std::sqrt(lambda) * a * u + squared_a * u) / (2 * lambda)};
}

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

  auto aa = (r1 * a * u + squared_a * u) / (std::complex<T>(2) * r1 * r1);
  auto aa2 = (r2 * a * u + squared_a * u) / (std::complex<T>(2) * r2 * r2);
  return {{r1, aa}, {r2, aa2}};
}

template<class T>
std::vector<std::pair<std::complex<T>,
                      Matrix<std::complex<T>>>> PowerMethodEigenvalues(
    const Matrix<T>& a) {
  if (!a.IsSquare()) {
    throw std::invalid_argument(
        "Matrix of size " + PairToString(a.Size()) + " is not square.");
  }
  auto squared_a = a * a;
  int n = a.Rows();

  auto m2_y = Matrix<T>(n, 1);
  m2_y(0) = 1;
  auto m2_u = m2_y / EuclideanNorm<T>(m2_y);
  auto m2_lambda = m2_u.ScalarProduct(a * m2_u);

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

  // while (EuclideanNorm<T>(a * u - lambda * u) > Matrix<T>::GetEps()) {
  for (int i = 0; i < 15; i++) {
    auto[p1, p2] = PowerMethodEigenvaluesComplexIteration(
        complex_a, complex_squared_a, complex_y, complex_u);
    auto r1 = p1.first;
    auto v1 = p1.second;
    auto r2 = p2.first;
    auto v2 = p2.second;
    std::cerr << r1 << ' ' << r2 << '\n';
    std::cerr << v1.Transposed() << v2.Transposed();
    auto
        error1 = norm(EuclideanNorm<std::complex<T>>(complex_a * v1 - r1 * v1));
    auto
        error2 = norm(EuclideanNorm<std::complex<T>>(complex_a * v2 - r2 * v2));
    std::cerr << error1 << ' ' << error2 << '\n';
    std::cerr << "---------\n";
  }
  // std::cerr << std::sqrt(m2_lambda) << '\n';
  // std::cerr << aa / EuclideanNorm<T>(aa);
  // std::cerr << aa2 / EuclideanNorm<T>(aa2);
  return {};
}

