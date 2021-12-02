#pragma once

#include "Matrix/matrix.h"
#include "Algebra/qr_decompose.h"

template<class T>
Matrix<T> MinimalSquareProblem(const Matrix<T>& a, const Matrix<T>& b) {
  auto[qr_a, qr_b] = QrDecompose(a, b);
  return SolveUxb(qr_a.SubMatrix(0, 0, a.Cols(), a.Cols()),
                  qr_b.SubMatrix(0, 0, a.Cols(), 1));
}
