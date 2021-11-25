#pragma once

template<class Matrix>
Matrix SolveLxb(const Matrix& l,
                const Matrix& b) {
  Matrix x(b.Rows(), 1);
  for (int i = 0; i < x.Rows(); i++) {
    auto sum = b.At(i, 0);
    for (int j = 0; j < i; j++) {
      sum -= x.At(j, 0) * l.At(i, j);
    }
    x.At(i, 0) = sum / l.At(i, i);
  }
  return x;
}

template<class Matrix>
Matrix SolveUxb(const Matrix& u,
                const Matrix& b) {
  Matrix x(b.Rows(), 1);
  for (int i = x.Rows() - 1; i >= 0; i--) {
    auto sum = b.At(i, 0);
    for (int j = x.Rows() - 1; j > i; j--) {
      sum -= x.At(j, 0) * u.At(i, j);
    }
    x.At(i, 0) = sum / u.At(i, i);
  }
  return x;
}

