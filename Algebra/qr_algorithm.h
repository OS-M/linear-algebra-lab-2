#pragma once

#include <complex>
#include "Matrix/matrix.h"
#include "rotations.h"
#include "eigenvalues.h"

template<class T>
int UnderDiagonalZeros(const Matrix<T>& a) {
  if (!a.IsSquare()) {
    throw std::invalid_argument(
        "Matrix of size " + PairToString(a.Size()) + " is not square.");
  }
  int ans = 0;
  for (int i = 0; i < a.Rows() - 1; i++) {
    if (std::abs(a(i + 1, i)) < Matrix<T>::GetEps()) {
      ans++;
    }
  }
  return ans;
}

template<class T>
bool DoDiagonalSquaresIntersect(const Matrix<T>& a) {
  if (!a.IsSquare()) {
    throw std::invalid_argument(
        "Matrix of size " + PairToString(a.Size()) + " is not square.");
  }
  for (int i = 0; i < a.Rows() - 2; i++) {
    if (std::abs(a(i + 1, i)) > Matrix<T>::GetEps() &&
        std::abs(a(i + 2, i + 1)) > Matrix<T>::GetEps()) {
      return true;
    }
  }
  return false;
}

template<class T>
std::vector<std::complex<T>> QrAlgorithm(Matrix<T> a,
                                         int* iters = nullptr,
                                         int max_iter = 1000) {
  if (!a.IsSquare()) {
    throw std::invalid_argument(
        "Matrix of size " + PairToString(a.Size()) + " is not square.");
  }
  int n = a.Rows();
  int iter;
  for (iter = 0; iter < max_iter; iter++) {
    std::vector<std::pair<T, T>> rotations;
    rotations.reserve(n);
    for (int i = 0; i < n - 1; i++) {
      auto[sin, cos] = GetRotationMatrix(a.SubMatrix(i, i, 2, 1));
      rotations.emplace_back(sin, cos);
      ApplyRotation(a, sin, cos, i);
    }

    for (int i = 0; i < n - 1; i++) {
      auto[sin, cos] = rotations[i];
      ApplyTransposedRotation(a, sin, cos, i);
    }

    if (DoDiagonalSquaresIntersect(a) ||
        UnderDiagonalZeros(a) < (n - 1) / 2) {
      continue;
    } else {
      break;
    }
  }
  if (iters) {
    *iters = iter + 1;
  }
  if (iter == max_iter) {
    if (iters) {
      *iters = -1;
    }
    return {};
  }

  std::vector<std::complex<T>> ans;
  for (int i = 0; i < n; i++) {
    if (i == n - 1 || std::abs(a(i + 1, i)) < Matrix<T>::GetEps()) {
      ans.emplace_back(a(i, i));
    } else {
      auto[e1, e2] = ExtractEigenvalues2x2(a.SubMatrix(i, i, 2, 2));
      ans.push_back(e1);
      ans.push_back(e2);
      i++;
    }
  }

  return ans;
}
