
#pragma once

#include <vector>
#include "Algebra/polynomial.h"

namespace __internal {

template<class T>
std::vector<T> FindRootsOfOddPolynomial(Polynomial<T> a, T eps) {
  return {};
}

}

template<class T>
std::vector<T> FindRoots(const Polynomial<T>& a, T eps) {
  if (a.size() == 2) {
    return {-a[1] / a[0]};
  }
  auto der = Derivative(a);
  std::vector<T> extremums{-1e18};
  for (auto it: FindRoots(der, eps)) {
    extremums.push_back(it);
  }
  extremums.push_back(1e18);
  std::vector<T> ans;
  for (int i = 1; i < extremums.size(); i++) {
    auto l = extremums[i - 1];
    auto r = extremums[i];
    bool rising = ValueIn(a, l) < ValueIn(a, r);
    while (!(std::abs(l - extremums[i]) < eps ||
        std::abs(r - extremums[i - 1]) < eps)) {
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
    auto val = ValueIn(a, (l + r) / 2);
    if (std::abs(ValueIn(a, (l + r) / 2)) < eps) {
      ans.push_back((l + r) / 2);
    }
  }
  return ans;
}
