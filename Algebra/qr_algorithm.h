#pragma once

#include "Matrix/matrix.h"
#include "rotations.h"

template <class T>
std::vector<T> QrAlgorithm(Matrix<T> a) {
  if (!a.IsSquare()) {
    throw std::invalid_argument(
        "Matrix of size " + PairToString(a.Size()) + " is not square.");
  }
  int n = a.Rows();
  for (int k = 0; k < 60; k++) {
    std::vector<std::pair<T, T>> rotations;
    for (int i = 0; i < n - 1; i++) {
      auto [sin, cos] = GetRotationMatrix(a.SubMatrix(i, i, 2, 1));
      rotations.emplace_back(sin, cos);
      ApplyRotation(a, sin, cos, i);
    }

    for (int i = 0; i < n - 1; i++) {
      auto [sin, cos] = rotations[i];
      ApplyTransposedRotation(a, sin, cos, i);
    }
  }
  std::cerr << a;

  return {};
}
