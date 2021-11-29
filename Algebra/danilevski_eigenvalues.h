#pragma once

#include "Matrix/matrix.h"

namespace __internal {

template<class T>
std::vector<T> DanilevskiPolynomial(const Matrix<T>& a) {
  std::vector<T> ans(a.Rows() + 1);
  ans[0] = 1;
  for (int i = 1; i < ans.size(); i++) {
    ans[i] = -a(0, i - 1);
  }
  if (a.Rows() % 2 == 1) {
    for (auto& it: ans) {
      it *= -1;
    }
  }
  return ans;
}

}

template<class T>
std::vector<std::vector<T>> DanilevskiPolynomial(Matrix<T>& a) {
  int last = 0;
  std::vector<std::vector<T>> ans;
  auto n = a.Rows();
  for (int i = 0; i < n - 1; i++) {
    if (std::abs(a(i + 1, i)) < Matrix<T>::GetEps()) {
      ans.push_back(__internal::DanilevskiPolynomial(
          a.SubMatrix(last, last, i - last + 1, i - last + 1)));
      last = i + 1;
    }
  }
  ans.push_back(__internal::DanilevskiPolynomial(
      a.SubMatrix(last, last, -1, -1)));
  return ans;
}
