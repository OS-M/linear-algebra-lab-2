#include <iostream>
#include "Matrix/matrix.h"
#include "Algebra/gauss.h"
#include "Algebra/euclidean_norm.h"
#include "Algebra/qr_hessenberg.h"

int main() {
  Matrix<double>::eps = 1e-6;

  auto a = DMatrix::Random(5, 5, 0, 1, 1);
  // auto b = DMatrix(3, 1, 1);
  // std::cout << a << b;
  std::cout << a.ToWolframString() << '\n';
  // auto x = GaussSolve(a, b).first;
  auto qr = QrHessenberg(a);
  // std::cout << qr;
  std::cout << qr.ToWolframString();
  // std::cout << QrHessenberg(a).ToWolframString();
  // std::cout << a.ToWolframString();
}
