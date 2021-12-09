
#pragma once

#include <vector>
#include <optional>
#include "Algebra/polynomial.h"

namespace __internal {

template<class T>
T RefineRootNewton(const Polynomial <T>& a,
                                  const Polynomial <T>& der,
                                  T x,
                                  T eps,
                                  int max_iters) {
  for (int i = 0; i < max_iters; ++i) {
    auto f_x = ValueIn(a, x);
    if (std::abs(f_x) < eps) {
      break;
    }
    auto der_x = ValueIn(der, x);
    x -= f_x / der_x;
  }
  return x;
}

template<class T>
std::optional<T> FindRootBinSearch(const Polynomial <T>& a, T l_, T r_, T eps) {
  auto l = l_;
  auto r = r_;
  bool rising = ValueIn(a, l) < ValueIn(a, r);
  int iters = 0;
  while (!(std::abs(l - r_) < eps || std::abs(r - l_) < eps)) {
    iters++;
    T mid = (l + r) / 2.;
    auto val = ValueIn(a, mid);
    if (std::abs(val) < eps) {
      break;
    }
    if (rising) {
      if (val < 0) {
        l = mid;
      } else {
        r = mid;
      }
    } else {
      if (val < 0) {
        r = mid;
      } else {
        l = mid;
      }
    }
  }
  // std::cerr << iters << '\n';
  if (std::abs(ValueIn(a, (l + r) / 2)) < eps) {
    return (l + r) / 2;
  }
  return {};
}

}

template<class T>
std::vector<T> FindRoots(Polynomial<T> a, T eps, T threshold) {
  auto der = Derivative(a);
  std::vector<T> ans;
  while (a.size() > 2) {
    auto x = __internal::RefineRootNewton(a, der, 1e9, eps, 1e4);
    auto root =
        __internal::FindRootBinSearch(a, x - threshold, x + threshold, eps);
    if (root.has_value()) {
      ans.push_back(root.value());
      a = DividePolynomial(a, {1, -root.value()});
      der = Derivative(a);
      // std::cerr << PolynomialToString(a);
    } else {
      break;
    }
  }
  if (a.size() == 2) {
    ans.emplace_back(-a[1] / a[0]);
  }
  return ans;
}
