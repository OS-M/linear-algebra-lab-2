#pragma once

#include "gauss_back_substitution.h"

template<class Matrix>
std::pair<Matrix, int> GaussSolve(Matrix a, Matrix b) {
  if (!a.IsSquare()) {
    throw std::invalid_argument(
        "Matrix of size " + PairToString(a.Size()) + " is not square.");
  }
  if (!b.IsColVector()) {
    throw std::invalid_argument(
        "Matrix of size " + PairToString(b.Size()) + " is not col.");
  }
  auto n = a.Rows();
  std::vector <size_t> reindex(n);
  std::iota(reindex.begin(), reindex.end(), 0);

  for (int i = 0; i < n; i++) {
    size_t index_of_max = i;
    for (int k = i + 1; k < n; k++) {
      if (a.At(reindex[index_of_max], i) < a.At(reindex[k], i)) {
        index_of_max = k;
      }
    }
    std::swap(reindex[i], reindex[index_of_max]);
    if (std::abs(a.At(reindex[i], i)) < Matrix::GetEps()) {
      continue;
    }
    for (int j = i + 1; j < n; j++) {
      if (std::abs(a.At(reindex[j], i)) < Matrix::GetEps()) {
        continue;
      }
      auto m = a.At(reindex[j], i) / a.At(reindex[i], i);
      a.SubMatrix(reindex[j], i, 1, -1) -=
          m * a.SubMatrix(reindex[i], i, 1, -1);
      b.At(reindex[j], 0) -= m * b.At(reindex[i], 0);
    }
  }

  auto swapped_a = a;
  auto swapped_b = b;
  int rank = 0;
  for (int i = 0; i < n; i++) {
    swapped_b.At(i, 0) = b.At(reindex[i], 0);
    for (int j = 0; j < n; j++) {
      swapped_a.At(i, j) = a.At(reindex[i], j);
    }
    if (std::fabs(swapped_a.At(i, i)) > Matrix::GetEps()) {
      rank++;
    }
  }
  return {SolveUxb(swapped_a, swapped_b), rank};
}
