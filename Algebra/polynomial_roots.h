
#pragma once

#include <vector>
#include "Algebra/polynomial.h"

namespace __internal {

template<class T>
std::vector<T> FindRootsOfOddPolynomial(Polynomial<T> a, T eps) {
  T l = -1e18;
  T r = 1e18;
  if (ValueIn(a, l) > ValueIn(a, r)) {
    for (auto& it: a) {
      it *= -1;
    }
  }
  while (std::abs(l - r) > eps) {
    T mid = (l + r) / 2.;
    auto val = ValueIn(a, mid);
    if (std::abs(val) < eps) {
      break;
    }
    if (val < 0) {
      l = mid;
    } else {
      r = mid;
    }
  }
  return {(l + r) / 2.};
}

}

template<class T>
std::vector<T> FindRoots(const Polynomial<T>& a, T eps) {
  if (a.size() % 2 == 0) {
    return __internal::FindRootsOfOddPolynomial(a, eps);
  }
  auto extremums = __internal::FindRootsOfOddPolynomial(Derivative(a), eps);
  return {};
}
