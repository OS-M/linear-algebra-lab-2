
#pragma once

#include <vector>
#include <optional>
#include "Algebra/polynomial.h"

namespace __internal {

template<class T>
std::optional<T> RefineRootNewton(const Polynomial <T>& a,
                                  const Polynomial <T>& der,
                                  T x,
                                  T eps,
                                  int max_iters) {
  for (int i = 0; i < max_iters; ++i) {
    auto f_x = ValueIn(a, x);
    if (std::abs(f_x) < eps) {
      // std::cerr << i << '\n';
      return x;
    }
    auto der_x = ValueIn(der, x);
    x -= f_x / der_x;
  }
  return {};
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

template<class T>
std::vector<T> FindRootsOfOddPolynomial(const Polynomial <T>& a, T eps) {
  if (a.size() == 2) {
    return {-a[1] / a[0]};
  }
  std::vector<T> ans;
  // for (int i = 1; i < extremums.size(); i++) {
  //   auto root = __internal::FindRootBinSearch(a,
  //                                             extremums[i - 1],
  //                                             extremums[i], eps);
  //   if (root.has_value()) {
  //     root = __internal::RefineRootNewton(a, der, root.value(), eps, 1e2);
  //     if (root.has_value()) {
  //       ans.push_back(root.value());
  //     }
  //   }
  // }
  return ans;
}

}

template<class T>
std::vector<T> FindRoots(const Polynomial <T>& a, T eps) {
  if (a.size() == 2) {
    return {-a[1] / a[0]};
  }
  auto der = Derivative(a);
  Normalize(der);
  std::vector<T> extremums{-1e4};
  for (auto it: FindRoots(der, eps)) {
    extremums.push_back(it);
  }
  extremums.push_back(1e4);
  std::vector<T> ans;
  for (int i = 1; i < extremums.size(); i++) {
    auto root = __internal::FindRootBinSearch(a,
                                              extremums[i - 1],
                                              extremums[i], eps);
    if (root.has_value()) {
      root = __internal::RefineRootNewton(a, der, root.value(), eps, 1e2);
      if (root.has_value()) {
        ans.push_back(root.value());
      }
    }
  }
  return ans;
}
