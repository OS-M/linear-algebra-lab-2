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
  for (int k = 0; k < 1; k++) {
    std::vector<Matrix<T>> rotations;
    for (int i = 0; i < n - 1; i++) {
      rotations.push_back(GetRotationMatrix(a.SubMatrix(i, i, 2, 1)));
      a.SubMatrix(i, i, 2, 2).CopyFrom(
          rotations.back() * a.SubMatrix(i, i, 2, 2));
      // a.SubMatrix(i, i + 1, 2, 1).CopyFrom(
      //     rotations.back() * a.SubMatrix(i, i + 1, 2, 1));
      // a(i, i) = std::sqrt(a(i, i) * a(i, i) + a(i + 1, i) * a(i + 1, i));
      // a(i + 1, i) = 0;
    }
    // std::cerr << a;
    for (int i = 0; i < n - 1; i++) {
      a.SubMatrix(i, i, 2, 2).CopyFrom(
          a.SubMatrix(i, i, 2, 2) * rotations[i].Transposed());
    }
    // break;
  }
  std::cerr << a;

  return {};
}
