#include <iostream>
#include "Matrix/matrix.h"
#include "Algebra/gauss.h"
#include "Algebra/euclidean_norm.h"
#include "Algebra/qr_hessenberg.h"
#include "Algebra/qr_algorithm.h"
#include "Algebra/power_method.h"

int main() {
  Matrix<double>::eps = 1e-6;

  {
    auto a = DMatrix::Random(3, 3, 1, 10, 1337);
    int k = 1;
    a(0, 0) = 5 * (k + 1);
    a(0, 1) = 4 * (k + 1);
    a(0, 2) = -2 * (1 + k);
    a(1, 0) = -6 - 5 * k;
    a(1, 1) = -5 - 4 * k;
    a(1, 2) = 2 * (1 + k);
    a(2, 0) = 2 * k;
    a(2, 1) = 2 * k;
    a(2, 2) = -1 - k;
    std::cout << a << a.ToWolframString();
    PowerMethodEigenvalues(a);
    DMatrix b(3, 1);
    b(0) = 1;
    b(1) = -1;
    b(2) = 1;
    std::cout << b / EuclideanNorm<double>(b);
  }
  return 0;
  {
    auto a = DMatrix::Random(5, 5, 1, 10, 1337);
    // auto b = DMatrix(3, 1, 1);
    // std::cout << a << b;
    // auto x = GaussSolve(a, b).first;
    // std::cout << x << a * x;

    std::cout << a << a.ToWolframString();
    auto qr = QrHessenberg(a);
    std::cout << qr << qr.ToWolframString();
    int iters = 0;
    auto values = QrAlgorithm(qr, &iters);
    for (auto val: values) {
      std::cout << val << '\n';
    }
    std::cout << iters;
    // std::cout << a.ToWolframString();
  }
}
