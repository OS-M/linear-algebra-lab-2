
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
  auto val_l = ValueIn(a, l);
  auto val_r = ValueIn(a, r);
  if ((val_l < 0 && val_r < 0) || (val_l > 0 && val_r > 0)) {
    std::cerr << "bad l r\n";
    return {};
  }
  bool rising = val_l > val_r;
  int iters = 0;
  while (!(std::abs(l - r_) < eps || std::abs(r - l_) < eps)
            && std::abs(l - r) > 1e-6) {
    iters++;
    T mid = (l + r) / 2.;
    auto val = ValueIn(a, mid);
    if (std::abs(val) < eps) {
      break;
    }
    if (rising) {
      if (val < 0) {
        r = mid;
      } else {
        l = mid;
      }
    } else {
      if (val < 0) {
        l = mid;
      } else {
        r = mid;
      }
    }
  }
  // std::cerr << iters << '\n';
  // if (std::abs(ValueIn(a, (l + r) / 2)) < eps) {
    return (l + r) / 2;
  // }
  return {};
}

}

template<class T>
std::vector<T> FindRoots(Polynomial<T> a, T eps, T threshold) {
  auto der = Derivative(a);
  std::vector<T> ans;
  while (a.size() > 2) {
    auto c = std::abs(*std::max_element(a.begin() + 1, a.end(), [] (T c1, T c2) {
      return std::abs(c1) < std::abs(c2);
    }));
    auto b = std::abs(*std::max_element(a.begin(), a.end() - 1, [] (T c1, T c2) {
      return std::abs(c1) < std::abs(c2);
    }));
    auto p = std::max(std::abs(a.back()) / (b + std::abs(a.back())),
                      1 + c / std::abs(a[0]));
    auto x = __internal::RefineRootNewton(a, der, p + 10, 100 * eps, 1e6);
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

template<class T>
std::vector<T> FindRoots2(const Polynomial <T>& a, T eps, int iter = 0) {
  if (a.size() == 2) {
    return {-a[1] / a[0]};
  }
  auto der = Derivative(a);
  for (auto& it: der) {
    it /= iter + 1;
  }
  std::cerr << PolynomialToString(a);
  std::vector<T> extremums{-1e18};
  for (auto it: FindRoots2(der, eps, iter + 1)) {
    extremums.push_back(it);
  }
  extremums.push_back(1e18);
  std::vector<T> ans;
  for (int i = 1; i < extremums.size(); i++) {
    auto root = __internal::FindRootBinSearch(a,
                                              extremums[i - 1],
                                              extremums[i], eps);
    if (root.has_value()) {
      // root = __internal::RefineRootNewton(a, der, root.value(), eps, 1e2);
      // if (root.has_value()) {
        ans.push_back(root.value());
      // }
    }
  }
  return ans;
}
