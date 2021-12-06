#include <iostream>
#include "Matrix/matrix.h"
#include "Algebra/gauss.h"
#include "Algebra/euclidean_norm.h"
#include "Algebra/reflections_hessenberg.h"
#include "Algebra/qr_algorithm.h"
#include "Algebra/power_iteration_method.h"
#include "Algebra/qr_decompose.h"
#include "Algebra/minimal_square_problem.h"
#include "Algebra/frobenius_form.h"
#include "Algebra/polynomial.h"
#include "Algebra/danilevski_eigenvalues.h"
#include "Algebra/polynomial_roots.h"
#include "TimeMeasurer/time_measurer.h"

DMatrix Matrix1() {
  return {{1, -2, 1, 0, -1, 1, -2, 2, 0, -2},
          {0, 2, 0, 0, 2, 1, -1, -1, -1, -2},
          {0, 1, 0, -1, 1, -1, 0, -1, 1, -1},
          {-2, -1, 2, -1, 0, 0, 0, 0, 1, 0},
          {1, -2, 0, 1, 0, -2, -1, 0, 2, 2},
          {-2, -2, 0, -2, 0, 1, 1, -2, 1, 1},
          {-1, -2, -1, -1, -2, -1, -2, 1, -1, 2},
          {-2, 1, 2, -2, 0, 2, 1, -1, -2, 2},
          {0, 1, 0, 1, 1, -2, 2, 0, 1, 1},
          {0, 0, 2, -1, -1, 0, -2, 2, -1, -1}};
}

DMatrix Matrix2() {
  return
      {{-1, 1, -1, 0, -1, 0, -1, 1, 1, -1, 0, -1, -1, 1, 0, 0, 1, 1, 1, 1},
       {-1, 0, -1, 1, -1, 0, 0, 0, 0, -1, 0, 0, -1, 1, 0, -1, 1, -1, -1, 0},
       {1, 0, -1, 1, 0, 1, -1, -1, -1, 0, -1, -1, 1, -1, 1, 1, -1, 1, -1, 0},
       {-1, 1, 0, 0, -1, 0, 0, -1, 0, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1, 0},
       {1, 0, -1, 0, 0, -1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, -1, 0, 0, 1},
       {0, 0, 0, 0, -1, 1, 1, 0, 0, 1, 1, 0, -1, 0, 1, 1, 0, 1, 0, 0},
       {-1, 0, 1, 1, 1, -1, -1, 0, -1, 1, -1, -1, -1, 0, -1, 0, 0, 0, -1, 1},
       {0, 0, -1, -1, 0, 1, 1, 1, 1, -1, 0, 0, -1, 1, 1, 1, 1, 0, 0, -1},
       {0, 0, 1, 1, 0, 1, 1, 0, 1, -1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1},
       {0, -1, 0, 0, 1, 0, -1, 0, -1, 0, -1, 0, -1, 0, 1, -1, 0, 0, 1, 1},
       {1, -1, 1, -1, -1, -1, 1, 0, -1, 0, 1, 1, -1, 0, 1, 1, 1, 0, 0, 0},
       {0, 1, 0, 0, -1, 0, 1, 0, 1, 0, 0, 1, 1, -1, -1, 0, -1, 1, 1, -1},
       {-1, -1, -1, -1, 0, 1, -1, 0, 0, -1, 0, 0, 0, 1, 1, 0, 0, 0, -1, 0},
       {-1, 0, 1, 0, -1, 0, 0, 1, -1, 1, 1, -1, 1, 1, 1, -1, 1, -1, -1, 0},
       {1, -1, 0, -1, -1, 0, -1, -1, 0, 0, 1, 0, 1, 1, -1, 1, 0, 0, -1, 0},
       {-1, -1, 1, 0, -1, 1, 1, -1, 1, 0, 0, -1, 1, -1, -1, 0, 0, 1, 1, 1},
       {0, 0, -1, 0, 0, 0, 0, -1, 1, 1, 0, -1, 1, -1, 0, 0, 0, -1, -1, 1},
       {-1, 0, -1, -1, -1, 1, 1, -1, 1, -1, 1, -1, 1, -1, 1, 1, 0, -1, 0, -1},
       {-1, 0, 1, 0, 0, 0, 0, -1, 1, -1, 1, -1, 0, -1, -1, 1, 0, 1, 0, 0},
       {0, -1, -1, 1, -1, 1, -1, -1, -1, 1, 1, -1, 0, -1, -1, 0, 1, 0, -1, -1}};
}

void Task1(double min, double max, int seed) {
  {
    int iters = 0;
    auto ans = PowerMethodEigenvalues(Matrix1(), &iters, 1000);
    std::cout << "Power iteration for matrix 1:\n";
    for (auto[e, v]: ans) {
      std::cout << e << '\n' << v;
    }
    std::cout << "Iters: " << iters << '\n';
    ans = PowerMethodEigenvalues(Matrix2(), &iters, 1000);
    std::cout << "Power iteration for matrix 2:\n";
    for (auto[e, v]: ans) {
      std::cout << e << '\n' << v;
    }
    std::cout << "Iters: " << iters << '\n';
  }

  std::vector<int> sizes{1000, 10000};
  std::vector<double> times;
  std::vector<int> iters;
  for (auto size: sizes) {
    times.push_back(0);
    iters.push_back(0);
    for (int t = 0; t < 5; t++) {
      auto a = DMatrix::Random(size, size, min, max, seed);
      TimeMeasurer time_measurer;
      int iter;
      PowerMethodEigenvalues(a, &iter, 1e2);
      seed ^= iter;
      times.back() += time_measurer.GetDuration();
      iters.back() += iter;
    }
    iters.back() /= 5;
    times.back() /= 5.;
    std::cout << times.back() << '\n';
    std::cout << iters.back() << '\n';
  }
}

template<class T>
void TestQrAlgorithm(const Matrix<T>& a) {
  std::cout << "Qr algorithm eigenvalues:\n";
  int iters = 0;
  auto qr_ans = QrAlgorithm(ReflectionsHessenberg(a), &iters, 12000);
  std::vector<std::complex<double>> complex_roots;
  complex_roots.reserve(qr_ans.size());
  for (auto r: qr_ans) {
    complex_roots.emplace_back(r);
  }
  auto vectors = FindEigenvectorsByValues(a.ToComplex(), complex_roots);
  for (int i = 0; i < vectors.size(); ++i) {
    std::cout << qr_ans[i] << '\n';
    if (vectors[i].Rows() != 0) {
      std::cout << vectors[i];
      std::cout << "Norm: " <<
                std::abs(EuclideanNorm<std::complex<T>>(
                    a.ToComplex() * vectors[i] - qr_ans[i] * vectors[i]));
    }
    std::cout << "\n#########\n";
  }
  std::cout << "Iters: " << iters << "\n===================\n\n\n";
}

template<class T>
void TestPowerMethod(const Matrix<T>& a) {
  int iters = 0;
  auto ans = PowerMethodEigenvalues(a, &iters, 1000, 20, 3);
  std::cout << "Power iteration eigenvalues and eigenvectors:\n";
  for (auto[e, v]: ans) {
    std::cout << e << '\n' << v;
    // std::cout << e << '\n';
    std::cout << "Norm: " <<
              std::abs(EuclideanNorm<std::complex<T>>(
                  a.ToComplex() * v - e * v));
    std::cout << "\n#########\n";
  }
  std::cout << "Iters: " << iters << "\n===================\n\n\n";
}

template<class T>
void TestDanilevskiMethod(const Matrix<T>& a) {
  auto polynomial = PolynomialMultiply(DanilevskiPolynomial(FrobeniusForm(a)));

  std::cout << "Danilevski Polynomial: " << PolynomialToString(polynomial);

  auto roots = FindRoots(polynomial, 1e-6);
  std::vector<std::complex<double>> complex_roots;
  complex_roots.reserve(roots.size());
  for (auto r: roots) {
    complex_roots.emplace_back(r);
  }

  auto vectors = FindEigenvectorsByValues(a.ToComplex(), complex_roots);
  for (int i = 0; i < vectors.size(); ++i) {
    if (vectors[i].Rows() == 0) {
      continue;
    }
    std::cout << roots[i] << '\n';
    std::cout << vectors[i];
    std::cout << "Norm: " <<
              std::abs(EuclideanNorm<std::complex<T>>(
                  a.ToComplex() * vectors[i] - roots[i] * vectors[i]));
    std::cout << "\n#########\n";
  }
  std::cout << "===================\n\n\n";
}

int main() {
  Matrix<double>::SetEps(1e-6, 6);
  Matrix<std::complex<double>>::SetEps(std::complex<double>(1e-6, 1e-6), 6);

  // {
  //   Polynomial<double> p{-1, -25, -59, -33};
  //   std::cout << PolynomialToString(p);
  //   for (auto r: FindRoots(p, 1e-6)) {
  //     std::cout << r << '\n';
  //   }
  //   return 0;
  // }
  // Task1(-100, 100, 228);

  {
    // auto a = DMatrix::RandomInts(3, 3, -5, 10, 228);
    auto a = DMatrix::RandomInts(3, 3, -5, 10, 1337);
    // auto a = DMatrix::Random(5, 5, -5, 10, 8541657);
    // auto a = DMatrix::Random(5, 5, -5, 10, 456468);
    // a = DMatrix{{1, 1, 1, 1, 1},
    //             {1, 1, 1, 1, 1},
    //             {0, 0, 1, 0, 0},
    //             {0, 0, 1, 0, 0},
    //             {0, 0, 0, 1, 0}};

    // a = DMatrix{{0, 1}, {1, 0}};

    int k = 2;
    a(0, 0) = 5 * (k + 1);
    a(0, 1) = 4 * (k + 1);
    a(0, 2) = -2 * (1 + k);
    a(1, 0) = -6 - 5 * k;
    a(1, 1) = -5 - 4 * k;
    a(1, 2) = 2 * (1 + k);
    a(2, 0) = 2 * k;
    a(2, 1) = 2 * k;
    a(2, 2) = -1 - k;

    // a = Matrix2();

    std::cout << a << "\n\n";

    TestQrAlgorithm(a);
    // TestDanilevskiMethod(a);
    TestPowerMethod(a);
  }
  return 0;

  {
    // auto a = DMatrix::RandomInts(3, 3, -5, 10, 228);
    // auto a = DMatrix::RandomInts(3, 3, -5, 10, 1337);
    auto a = DMatrix::Random(5, 5, -5, 10, 8541657);
    // a = DMatrix{{1, 1, 1, 1, 1},
    //             {1, 1, 1, 1, 1},
    //             {0, 0, 1, 0, 0},
    //             {0, 0, 1, 0, 0},
    //             {0, 0, 0, 1, 0}};

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
    auto complex_a = a.ToComplex();
    std::cout << FrobeniusForm(a).ToWolframString();
    auto p = DanilevskiPolynomial(FrobeniusForm(a));
    std::vector<double> full_p(1, 1);
    for (auto p_: p) {
      // std::cout << PolynomialToString(p_);
      full_p = PolynomialMultiply(full_p, p_);
    }
    // Normalize(full_p);
    std::cout << PolynomialToString(full_p);
    std::vector<std::complex<double>> complex;
    auto roots = FindRoots(full_p, 1e-6);
    for (auto r: roots) {
      complex.emplace_back(r);
    }
    auto vectors = FindEigenvectorsByValues(complex_a, complex);
    for (int i = 0; i < vectors.size(); ++i) {
      std::cout << roots[i] << '\n';
      std::cout << vectors[i];
      std::cout << EuclideanNorm<std::complex<double>>(
          complex_a * vectors[i] - roots[i] * vectors[i]) << '\n';
    }

    int iters = 0;
    auto ans = PowerMethodEigenvalues(a, &iters, 1000);
    std::cout << "Power iteration:\n";
    for (auto[e, v]: ans) {
      std::cout << e << '\n' << v;
    }
    std::cout << "Iters: " << iters << '\n';
    std::cout << "Qr:\n";
    auto qr_ans = QrAlgorithm(ReflectionsHessenberg(a), &iters, 12000);
    for (auto val: qr_ans) {
      std::cout << val << '\n';
    }
    std::cout << "Iters: " << iters << '\n';
  }
  // std::vector<double> p{1, 2, 3};
  // std::vector<double> p2{1, 2};
  // std::cout << PolynomialToString(p) << PolynomialToString(p2)
  //           << PolynomialToString(PolynomialMultiply(p, p2));
  return 0;

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
    auto qr_ans = QrAlgorithm(ReflectionsHessenberg(a), &iters, 12000);
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
    auto qr = ReflectionsHessenberg(a);
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
