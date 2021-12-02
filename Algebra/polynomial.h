#pragma once

#include <vector>

template<class T>
using Polynomial = std::vector<T>;

template<class T>
Polynomial<T> PolynomialMultiply(const Polynomial<T>& a,
                                 const Polynomial<T>& b) {
  Polynomial<T> c(a.size() + b.size() - 1);
  for (int i = 0; i < a.size(); i++) {
    for (int j = 0; j < b.size(); j++) {
      c[i + j] += a[i] * b[j];
    }
  }
  return c;
}

template<class T>
std::string PolynomialToString(const Polynomial<T>& a,
                               const std::string& x = "x") {
  std::stringstream ss;
  for (int i = 0; i < a.size(); i++) {
    ss << a[i] << " * " << x << "^" << a.size() - i - 1 << " + ";
  }
  auto ans = ss.str();
  return ans.substr(0, ans.size() - 3) + "\n";
}

template<class T>
Polynomial<T> Derivative(const Polynomial<T>& a) {
  Polynomial<T> ans = a;
  ans.pop_back();
  for (int i = 0; i < ans.size(); i++) {
    ans[i] *= (ans.size() - i);
  }
  return ans;
}

template<class T>
T ValueIn(const Polynomial<T>& a, T x) {
  T ans = 0;
  T x_pow = 1;
  for (int i = a.size() - 1; i >= 0; i--) {
    ans += x_pow * a[i];
    x_pow *= x;
  }
  return ans;
}
