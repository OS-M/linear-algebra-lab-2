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
  // while (EuclideanNorm<T>(a * u - lambda * u) > Matrix<T>::eps) {
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
  // while (EuclideanNorm<T>(a * u - lambda * u) > Matrix<T>::eps) {
  for (int i = 0; i < 15; i++) {
    y.Assign(a2 * u);
    u.Assign(y / EuclideanNorm<T>(y));
    lambda = (u.ScalarProduct(a2 * u));
  }
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
  // while (EuclideanNorm<T>(a * u - lambda * u) > Matrix<T>::eps) {
  for (int g = 0; g < 150; g++) {
    y.Assign(a * u);
    u.Assign(y / EuclideanNorm<T>(y));
    lambda = (u.ScalarProduct(a * u));

    Matrix<T> l(u.Rows(), 2);
    l.Col(0).Assign(u);
    l.Col(1).Assign(a * u);
    auto r = -1 * a * a * u;
    auto c = MinimalSquareProblem(l, r);
    auto [r1, r2] = SolveQuadraticEquation<T>(1, c(1), c(0));
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
    auto aa = (r1 * a_c * u_c + a_c * a_c * u_c) / (std::complex<T>(2) * r1 * r1);
    // auto aa2 = (lambda * a * u + a * a * u) / (2 * lambda * lambda);
    std::cerr << norm(EuclideanNorm<std::complex<T>>(a_c * aa - r1 * aa)) << '\n';
    // std::cerr << aa2 / EuclideanNorm<T>(aa2);
    std::cerr << "--------\n";
  }
  // std::cerr << lambda << '\n' << u;
  // Matrix<T> l(u.Rows(), 2);
  // l.Col(0).Assign(u);
  // l.Col(1).Assign(a * u);
  // auto r = a * a * u;
  // auto c = MinimalSquareProblem(l, r);
  // std::cerr << c;
  // auto [r1, r2] = SolveQuadraticEquation<T>(1, c(1), c(0));
  // std::cerr << r1 << ' ' << r2 << '\n';
  return {};
}
