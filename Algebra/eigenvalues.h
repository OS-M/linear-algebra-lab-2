#pragma once

#include "Matrix/matrix.h"

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
