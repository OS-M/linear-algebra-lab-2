#pragma once

#include "Matrix/matrix.h"

enum class RowOperation {
  kSwap,
  kAdd,
  kMultiply,
};

template<class T>
Matrix<T> FrobeniusForm(Matrix<T> a,
                        std::vector<std::vector<std::tuple
                            <RowOperation, int, int, T>>>* operations = nullptr) {
  if (!a.IsSquare()) {
    throw std::invalid_argument(
        "Matrix of size " + PairToString(a.Size()) + " is not square.");
  }
  auto n = a.Rows();

  std::vector<std::tuple<RowOperation, int, int, T>> now_opers;

  for (int i = n - 1; i > 0; i--) {
    int row = i;
    int col = i - 1;
    int index_of_max = 0;

    for (int j = 0; j < i; j++) {
      if (std::abs(a(row, index_of_max)) < std::abs(a(row, j))) {
        index_of_max = j;
      }
    }
    if (i - 1 != index_of_max) {
      Matrix<T> r;
      r.CopyFrom(a.Col(i - 1));
      a.Col(i - 1).Assign(a.Col(index_of_max));
      a.Col(index_of_max).Assign(r);

      r.CopyFrom(a.Row(i - 1));
      a.Row(i - 1).Assign(a.Row(index_of_max));
      a.Row(index_of_max).Assign(r);

      now_opers.emplace_back(RowOperation::kSwap, i - 1, index_of_max, 0);
    }

    if (std::abs(a(row, col)) < Matrix<T>::GetEps()) {
      std::cerr << std::fixed << std::setprecision(30) << a(row, col) << '\n';
      if (operations) {
        operations->push_back(now_opers);
        now_opers.clear();
      }
      continue;
    }
    auto d = a(row, col);
    now_opers.emplace_back(RowOperation::kMultiply, col, col, 1. / d);
    a.Col(col).SubMatrix(0, 0, i + 1, -1) /= d;
    a.Row(col) *= d;

    for (int j = 0; j < n; j++) {
      if (col == j) {
        continue;
      }
      auto d = a(row, j);
      now_opers.emplace_back(RowOperation::kAdd, col, j, -d);
      for (int k = 0; k < i + 1; k++) {
        a(k, j) -= a(k, col) * d;
      }
      for (int k = 0; k < n; k++) {
        a(col, k) += a(j, k) * d;
      }
    }
  }
  operations->push_back(now_opers);
  return a;
}

template<class T>
std::vector<Matrix<T>> EigenVectorsForFrobeniusForm(
    int n,
    int shift,
    int cur_size,
    const std::vector<std::tuple<RowOperation, int, int, T>>& operations_,
    const std::vector<T>& eigenvalues) {
  std::vector<Matrix<T>> ans;
  auto operations = operations_;
  std::cerr << operations_.size() << '\n';
  std::reverse(operations.begin(), operations.end());
  for (auto eigenvalue: eigenvalues) {
    Matrix<T> v(n, 1);
    v(shift + cur_size - 1) = 1;
    for (int i = shift + cur_size - 2; i >= shift; i--) {
      v(i) = v(i + 1) * eigenvalue;
    }
    for (auto[operation, k1, k2, k3] : operations) {
      switch (operation) {
        case RowOperation::kSwap:
          std::swap(v(k1), v(k2));
          break;

        case RowOperation::kAdd:
          v(k1) += k3 * v(k2);
          break;

        case RowOperation::kMultiply:
          v(k1) *= k3;
          break;
      }
    }
    ans.push_back(v);
  }
  return ans;
}
