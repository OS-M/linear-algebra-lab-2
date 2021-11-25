#include <iostream>
#include "Matrix/matrix.h"

int main() {
  Matrix<double>::eps = 1e-6;

  auto a = DMatrix::Ones(3);
  a(1, 2) = 1;
  auto b = a * a.Transposed() * 3;
  std::cout << a << b;
  b.Col(1) *= 2;
  b.Row(0) *= 2;
  auto c = b.SubMatrix(2, 1, 1, 2).Transposed();
}
