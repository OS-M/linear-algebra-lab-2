#include <iostream>
#include "Matrix/matrix.h"
#include "Algebra/gauss.h"
#include "Algebra/euclidean_norm.h"
#include "Algebra/qr_hessenberg.h"
#include "Algebra/qr_algorithm.h"

int main() {
  Matrix<double>::eps = 1e-6;

  auto a = DMatrix::Random(5, 5, 1, 10, 5);
  // auto b = DMatrix(3, 1, 1);
  // std::cout << a << b;
  // auto x = GaussSolve(a, b).first;
  // std::cout << x << a * x;

  std::cout << a << a.ToWolframString();
  auto qr = QrHessenberg(a);
  std::cout << qr << qr.ToWolframString();
  auto values = QrAlgorithm(qr);
  for (auto val: values) {
    std::cout << val << '\n';
  }
  // std::cout << a.ToWolframString();
}
