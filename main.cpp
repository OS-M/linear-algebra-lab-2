#include <iostream>
#include <fstream>
#include <mutex>
#include <shared_mutex>
#include <thread>
#include "Matrix/matrix.h"
#include "Algebra/gauss.h"
#include "Algebra/euclidean_norm.h"
#include "Algebra/hessenberg_form.h"
#include "Algebra/qr_algorithm.h"
#include "Algebra/power_iteration_method.h"
#include "Algebra/qr_decompose.h"
#include "Algebra/minimal_square_problem.h"
#include "Algebra/frobenius_form.h"
#include "Algebra/polynomial.h"
#include "Algebra/danilevski_eigenvalues.h"
#include "Algebra/polynomial_roots.h"
#include "Plot/plot.h"
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
  // {
  //   int iters = 0;
  //   auto ans = PowerMethodEigenvalues(Matrix1(), &iters, 1000);
  //   std::cout << "Power iteration for matrix 1:\n";
  //   for (auto[e, v]: ans) {
  //     std::cout << e << '\n' << v;
  //   }
  //   std::cout << "Iters: " << iters << '\n';
  //   ans = PowerMethodEigenvalues(Matrix2(), &iters, 1000);
  //   std::cout << "Power iteration for matrix 2:\n";
  //   for (auto[e, v]: ans) {
  //     std::cout << e << '\n' << v;
  //   }
  //   std::cout << "Iters: " << iters << '\n';
  // }

  std::vector<int> sizes{10, 50, 100, 1000};
  int tests_count = 10;
  int max_iter = 1e4;
  std::vector<std::vector<DMatrix>> v;
  std::shared_mutex mutex;
  int thread_num = 12;
  std::cout << "Generating matrices...\n";
  auto thread_main = [&](int size, int id) {
    while (true) {
      auto a = DMatrix::Random(size, size, min, max, seed + id);
      int iter = -1;
      PowerMethodEigenvalues(a, &iter, max_iter, 20, 5, 0);
      mutex.lock_shared();
      if (v.back().size() >= tests_count) {
        mutex.unlock_shared();
        break;
      }
      mutex.unlock_shared();
      if (iter > 0) {
        mutex.lock();
        v.back().push_back(a);
        std::cout << "!\n";
        mutex.unlock();
      }
    }
  };
  for (auto size: sizes) {
    v.emplace_back();
    std::vector<std::thread> threads;
    threads.reserve(thread_num);
    std::cout << size << '\n';
    for (int i = 0; i < thread_num; i++) {
      threads.emplace_back(thread_main, size, i);
    }
    for (auto& thread: threads) {
      thread.join();
    }
    while (v.back().size() > tests_count) {
      v.back().pop_back();
    }
    std::cout << v.back().size() << '\n';
  }

  Plot times_plot("Times", "size", "time", sizes);
  Plot iters_plot("Iters", "size", "iters", sizes);

  std::vector<std::string> methods{"Mod 1", "Mod 2", "Mod 3", "Auto"};
  for (int method_ind = 0; method_ind < methods.size(); ++method_ind) {
    PlotLine times_line(methods[method_ind]);
    PlotLine iters_line(methods[method_ind]);
    for (int i = 0; i < sizes.size(); i++) {
      double time = 0;
      int iters = 0;
      for (const auto& m: v[i]) {
        TimeMeasurer time_measurer;
        int iter;
        PowerMethodEigenvalues(m, &iter, 2 * max_iter, 20, 5, method_ind);
        if (iter < 0) {
          std::cerr << "Not converge" << '\n';
          iter = 1000 * max_iter;
        }
        iters += iter;
        time += time_measurer.GetDuration();
      }
      iters /= tests_count;
      time /= tests_count;
      times_line.AddValue(sizes[i], time);
      iters_line.AddValue(sizes[i], iters);
      std::cout << sizes[i] << '\n';
    }
    times_plot.AddPlotLine(times_line);
    iters_plot.AddPlotLine(iters_line);
  }
  std::ofstream out("../task1_plot.txt");
  out << times_plot.ToString() << iters_plot.ToString();
}

void Task1_(double min, double max, int seed) {
  std::vector<int> sizes{10, 50, 100, 200, 500};
  int tests_count = 10;
  int max_iter = 1e4;
  std::vector<std::vector<DMatrix>> v;
  std::shared_mutex mutex;
  int thread_num = 12;
  std::cout << "Generating matrices...\n";
  auto thread_main = [&](int size, int id) {
    while (true) {
      auto a = DMatrix::Random(size, size, min, max, seed + id);
      int iter = -1;
      PowerMethodEigenvalues(a, &iter, max_iter, 20, 5, 1);
      mutex.lock_shared();
      if (v.back().size() >= tests_count) {
        mutex.unlock_shared();
        break;
      }
      mutex.unlock_shared();
      if (iter > 0) {
        mutex.lock();
        v.back().push_back(a);
        std::cout << "!\n";
        mutex.unlock();
      }
    }
  };
  for (auto size: sizes) {
    v.emplace_back();
    std::vector<std::thread> threads;
    threads.reserve(thread_num);
    std::cout << size << '\n';
    for (int i = 0; i < thread_num; i++) {
      threads.emplace_back(thread_main, size, i);
    }
    for (auto& thread: threads) {
      thread.join();
    }
    while (v.back().size() > tests_count) {
      v.back().pop_back();
    }
    std::cout << v.back().size() << '\n';
  }

  Plot times_plot("Times", "size", "time", sizes);
  Plot iters_plot("Iters", "size", "iters", sizes);

  std::vector<std::string> methods{"Mod 1", "Mod 2", "Mod 3", "Auto"};
  for (int method_ind = 1; method_ind < methods.size(); ++method_ind) {
    PlotLine times_line(methods[method_ind]);
    PlotLine iters_line(methods[method_ind]);
    for (int i = 0; i < sizes.size(); i++) {
      double time = 0;
      int iters = 0;
      for (const auto& m: v[i]) {
        TimeMeasurer time_measurer;
        int iter;
        PowerMethodEigenvalues(m, &iter, 5 * max_iter, 20, 5, method_ind);
        if (iter < 0) {
          std::cerr << "Not converge " << methods[method_ind] << '\n';
          iter = 1000 * max_iter;
        }
        iters += iter;
        time += time_measurer.GetDuration();
      }
      iters /= tests_count;
      time /= tests_count;
      times_line.AddValue(sizes[i], time);
      iters_line.AddValue(sizes[i], iters);
      std::cout << sizes[i] << '\n';
    }
    times_plot.AddPlotLine(times_line);
    iters_plot.AddPlotLine(iters_line);
  }
  std::ofstream out("../task1_plot2.txt");
  out << times_plot.ToString() << iters_plot.ToString();
}

void Task1__(double min, double max, int seed) {
  int size = 50;
  int count = 1000;
  int max_iter = 2000;
  int thread_num = 12;
  std::vector<std::thread> threads;
  std::vector<std::vector<int>> thread_ans;
  auto thread_main = [&](
      std::vector<int>& v, int id, int algorithm_id) {
    for (int i = 0; i < count; i++) {
      auto a = DMatrix::Random(size, size, min, max, seed + id);
      int iter = 0;
      PowerMethodEigenvalues(a, &iter, max_iter, 10, 5, algorithm_id);
      if (iter < 0) {
        iter = max_iter;
      }
      iter = std::min(iter, max_iter);
      v[iter]++;
      if ((i + 1) % 100 == 0) {
        std::cout << id << ": " << i << '\n';
      }
    }
  };
  thread_ans.resize(thread_num, std::vector<int>(max_iter + 1));
  for (int i = 0; i < thread_num; i++) {
    threads.emplace_back(thread_main, std::ref(thread_ans[i]), i + 1, 0);
  }
  for (auto& thread: threads) {
    thread.join();
  }
  std::vector<int> ans(max_iter + 1);
  for (int i = 0; i < ans.size(); i++) {
    for (auto& thread_answer: thread_ans) {
      ans[i] += thread_answer[i];
    }
  }
  std::vector<int> xs(max_iter + 1);
  std::iota(xs.begin(), xs.end(), 0);
  Plot plot("Times", "size", "time", xs);
  PlotLine line("Auto");
  for (int i = 0; i < ans.size(); i++) {
    line.AddValue(i, ans[i]);
  }
  plot.AddPlotLine(line);
  std::ofstream out("../task1_plot3.txt");
  out << plot.ToString();
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
      // std::cout << vectors[i];
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

  auto roots = FindRoots(polynomial, 1e-6, 1.);
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
    // std::cout << vectors[i];
    std::cout << "Norm: " <<
              std::abs(EuclideanNorm<std::complex<T>>(
                  a.ToComplex() * vectors[i] - roots[i] * vectors[i]));
    std::cout << "\n#########\n";
  }
  std::cout << "===================\n\n\n";
}

int main() {
  auto eps = 1e-6;
  auto prec = 6;
  Matrix<double>::SetEps(eps, prec);
  Matrix<std::complex<double>>::SetEps(std::complex<double>(eps, eps), prec);

  // {
  //   Polynomial<double> a{1, 226.99, -310.27, -17936.9, 68579.4, -684};
  //   std::cout << PolynomialToString(a);
  //   for (auto root: FindRoots(a, eps)) {
  //     std::cout << root << '\n';
  //   }
  //   return 0;
  // }
  // Task1(-1e6, 1e6, 8917293);
  // Task1_(-1e6, 1e6, 8917293);
  // Task1__(-1e6, 1e6, 8917293);
  // return 0;

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
    a = Matrix1();

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

    // a = Matrix2();

    std::cout << a << "\n\n";

    TestQrAlgorithm(a);
    TestDanilevskiMethod(a);
    TestPowerMethod(a);
  }
}
