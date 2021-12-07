#pragma once

template<class Matrix>
Matrix SolveLxb(const Matrix& l,
                const Matrix& b) {
  if (!l.IsSquare()) {
    throw std::invalid_argument(
        "Matrix of size " + PairToString(l.Size()) + " is not square.");
  }
  if (!b.IsColVector()) {
    throw std::invalid_argument(
        "Matrix of size " + PairToString(b.Size()) + " is not col.");
  }
  if (b.Rows() != l.Rows()) {
    throw std::invalid_argument(
        "Bad matrix sizes: " + PairToString(l.Size()) + " and "
            + PairToString(b.Size()));
  }
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
  if (!u.IsSquare()) {
    throw std::invalid_argument(
        "Matrix of size " + PairToString(u.Size()) + " is not square.");
  }
  if (!b.IsColVector()) {
    throw std::invalid_argument(
        "Matrix of size " + PairToString(b.Size()) + " is not col.");
  }
  if (b.Rows() != u.Rows()) {
    throw std::invalid_argument(
        "Bad matrix sizes: " + PairToString(u.Size()) + " and "
            + PairToString(b.Size()));
  }
  Matrix x(b.Rows(), 1);
  for (int i = x.Rows() - 1; i >= 0; i--) {
    auto sum = b.At(i, 0);
    for (int j = x.Rows() - 1; j > i; j--) {
      sum -= x.At(j, 0) * u.At(i, j);
    }
    if (std::abs(u.At(i, i)) >= 5 * std::abs(Matrix::GetEps())) {
      x.At(i, 0) = sum / u.At(i, i);
    } else {
      x.At(i, 0) = 1;
    }
  }
  return x;
}

