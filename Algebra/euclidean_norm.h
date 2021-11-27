#pragma once

#include <complex>

template<class T, class Matrix>
T EuclideanNorm(const Matrix& a) {
  T res = T();
  for (int i = 0; i < a.Rows(); i++) {
    for (int j = 0; j < a.Cols(); j++) {
      res += a.At(i, j) * a.At(i, j);
    }
  }
  return std::sqrt(res);
}
