#include <iostream>
#include "Matrix/matrix.h"
#include "Algebra/gauss.h"
#include "Algebra/euclidean_norm.h"
#include "Algebra/qr_hessenberg.h"
#include "Algebra/qr_algorithm.h"
#include "Algebra/power_iteration_method.h"
#include "Algebra/qr_decompose.h"
#include "Algebra/minimal_square_problem.h"
#include "Algebra/frobenius_form.h"
#include "Algebra/polynomial.h"
#include "Algebra/danilevski_eigenvalues.h"

int main() {
  Matrix<double>::SetEps(1e-6, 6);
  Matrix<std::complex<double>>::SetEps(std::complex<double>(1e-6, 1e-6), 6);

  {
    // auto a = DMatrix::RandomInts(3, 3, -5, 10, 228);
    auto a = DMatrix::RandomInts(3, 3, -5, 10, 1337);
    // a = DMatrix{{1, 1, 1, 1, 1},
    //             {1, 1, 1, 1, 1},
    //             {0, 0, 1, 0, 0},
    //             {0, 0, 1, 0, 0},
    //             {0, 0, 0, 1, 0}};
    // a = DMatrix{{1, -2, 1, 0, -1, 1, -2, 2, 0, -2},
    //             {0, 2, 0, 0, 2, 1, -1, -1, -1, -2},
    //             {0, 1, 0, -1, 1, -1, 0, -1, 1, -1},
    //             {-2, -1, 2, -1, 0, 0, 0, 0, 1, 0},
    //             {1, -2, 0, 1, 0, -2, -1, 0, 2, 2},
    //             {-2, -2, 0, -2, 0, 1, 1, -2, 1, 1},
    //             {-1, -2, -1, -1, -2, -1, -2, 1, -1, 2},
    //             {-2, 1, 2, -2, 0, 2, 1, -1, -2, 2},
    //             {0, 1, 0, 1, 1, -2, 2, 0, 1, 1},
    //             {0, 0, 2, -1, -1, 0, -2, 2, -1, -1}};
    int k = 2;
    // a(0, 0) = 5 * (k + 1);
    // a(0, 1) = 4 * (k + 1);
    // a(0, 2) = -2 * (1 + k);
    // a(1, 0) = -6 - 5 * k;
    // a(1, 1) = -5 - 4 * k;
    // a(1, 2) = 2 * (1 + k);
    // a(2, 0) = 2 * k;
    // a(2, 1) = 2 * k;
    // a(2, 2) = -1 - k;
    std::cout << a << a.ToWolframString();
    auto form = FrobeniusForm(a);
    std::cout << form << form.ToWolframString();
    auto p = DanilevskiPolynomial(form);
    std::vector<double> full_p(1, 1);
    for (auto p_: p) {
      std::cout << PolynomialToString(p_);
      full_p = PolynomialMultiply(full_p, p_);
    }
    std::cout << PolynomialToString(full_p);
  }
  // std::vector<double> p{1, 2, 3};
  // std::vector<double> p2{1, 2};
  // std::cout << PolynomialToString(p) << PolynomialToString(p2)
  //           << PolynomialToString(PolynomialMultiply(p, p2));
  // return 0;

  {
    auto a = DMatrix::RandomInts(3, 3, -5, 10, 1337);
    // a = DMatrix{{1, -2, 1, 0, -1, 1, -2, 2, 0, -2},
    //             {0, 2, 0, 0, 2, 1, -1, -1, -1, -2},
    //             {0, 1, 0, -1, 1, -1, 0, -1, 1, -1},
    //             {-2, -1, 2, -1, 0, 0, 0, 0, 1, 0},
    //             {1, -2, 0, 1, 0, -2, -1, 0, 2, 2},
    //             {-2, -2, 0, -2, 0, 1, 1, -2, 1, 1},
    //             {-1, -2, -1, -1, -2, -1, -2, 1, -1, 2},
    //             {-2, 1, 2, -2, 0, 2, 1, -1, -2, 2},
    //             {0, 1, 0, 1, 1, -2, 2, 0, 1, 1},
    //             {0, 0, 2, -1, -1, 0, -2, 2, -1, -1}};
    // a = DMatrix{
    //     {-1, 1, -1, 0, -1, 0, -1, 1, 1, -1, 0, -1, -1, 1, 0, 0, 1, 1, 1, 1},
    //     {-1, 0, -1, 1, -1, 0, 0, 0, 0, -1, 0, 0, -1, 1, 0, -1, 1, -1, -1, 0},
    //     {1, 0, -1, 1, 0, 1, -1, -1, -1, 0, -1, -1, 1, -1, 1, 1, -1, 1, -1, 0},
    //     {-1, 1, 0, 0, -1, 0, 0, -1, 0, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 0},
    //     {1, 0, -1, 0, 0, -1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, -1, 0, 0, 1},
    //     {0, 0, 0, 0, -1, 1, 1, 0, 0, 1, 1, 0, -1, 0, 1, 1, 0, 1, 0, 0},
    //     {-1, 0, 1, 1, 1, -1, -1, 0, -1, 1, -1, -1, -1, 0, -1, 0, 0, 0, -1, 1},
    //     {0, 0, -1, -1, 0, 1, 1, 1, 1, -1, 0, 0, -1, 1, 1, 1, 1, 0, 0, -1},
    //     {0, 0, 1, 1, 0, 1, 1, 0, 1, -1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1},
    //     {0, -1, 0, 0, 1, 0, -1, 0, -1, 0, -1, 0, -1, 0, 1, -1, 0, 0, 1, 1},
    //     {1, -1, 1, -1, -1, -1, 1, 0, -1, 0, 1, 1, -1, 0, 1, 1, 1, 0, 0, 0},
    //     {0, 1, 0, 0, -1, 0, 1, 0, 1, 0, 0, 1, 1, -1, -1, 0, -1, 1, 1, -1},
    //     {-1, -1, -1, -1, 0, 1, -1, 0, 0, -1, 0, 0, 0, 1, 1, 0, 0, 0, -1, 0},
    //     {-1, 0, 1, 0, -1, 0, 0, 1, -1, 1, 1, -1, 1, 1, 1, -1, 1, -1, -1, 0},
    //     {1, -1, 0, -1, -1, 0, -1, -1, 0, 0, 1, 0, 1, 1, -1, 1, 0, 0, -1, 0},
    //     {-1, -1, 1, 0, -1, 1, 1, -1, 1, 0, 0, -1, 1, -1, -1, 0, 0, 1, 1, 1},
    //     {0, 0, -1, 0, 0, 0, 0, -1, 1, 1, 0, -1, 1, -1, 0, 0, 0, -1, -1, 1},
    //     {-1, 0, -1, -1, -1, 1, 1, -1, 1, -1, 1, -1, 1, -1, 1, 1, 0, -1, 0, -1},
    //     {-1, 0, 1, 0, 0, 0, 0, -1, 1, -1, 1, -1, 0, -1, -1, 1, 0, 1, 0, 0},
    //     {0, -1, -1, 1, -1, 1, -1, -1, -1, 1, 1, -1, 0, -1, -1, 0, 1, 0, -1,
    //      -1}};
    int k = 2;
    // a(0, 0) = 5 * (k + 1);
    // a(0, 1) = 4 * (k + 1);
    // a(0, 2) = -2 * (1 + k);
    // a(1, 0) = -6 - 5 * k;
    // a(1, 1) = -5 - 4 * k;
    // a(1, 2) = 2 * (1 + k);
    // a(2, 0) = 2 * k;
    // a(2, 1) = 2 * k;
    // a(2, 2) = -1 - k;

    std::cout << a << a.ToWolframString();

    int iters = 0;
    auto ans = PowerMethodEigenvalues(a, &iters, 1000);
    std::cout << "Power iteration:\n";
    for (auto[e, v]: ans) {
      std::cout << e << '\n' << v;
    }
    std::cout << "Iters: " << iters << '\n';
    std::cout << "Qr:\n";
    auto qr_ans = QrAlgorithm(QrHessenberg(a), &iters, 12000);
    for (auto val: qr_ans) {
      std::cout << val << '\n';
    }
    std::cout << "Iters: " << iters << '\n';
    // DMatrix b{{1.36, 0.88, 1}};
    // std::cout << b.Transposed() / EuclideanNorm<double>(b);
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

  // {
  //   DMatrix a{{1, 2},
  //             {2, 3},
  //             {3, 4},
  //             {10, 100}};
  //   DMatrix b{{1, 1, 3, 100}};
  //   std::cout << a.Transposed().ToWolframString() << b.ToWolframString();
  //   b = b.Transposed();
  //   // auto gauss_x = GaussSolve(a, b).first;
  //   // std::cout << a << b << gauss_x;
  //   // std::cout << a * gauss_x - b;
  //   auto solve = MinimalSquareProblem(a, b);
  //   auto gauss_solve = GaussSolve(a.Transposed() * a, a.Transposed() * b).first;
  //   std::cout << solve << gauss_solve;
  //   std::cout << a * solve - b;
  //   std::cout << a * gauss_solve - b;
  // }
  // return 0;
}
