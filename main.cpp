#include <iostream>
#include "Matrix/matrix.h"
#include "Algebra/gauss.h"

int main() {
  Matrix<double>::eps = 1e-6;

  auto a = DMatrix::Random(3, 3, 0, 1);
  auto b = DMatrix(3, 1, 1);
  std::cout << a << b;
  auto x = GaussSolve(a, b).first;
  std::cout << x << a * x;
}
