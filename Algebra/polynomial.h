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
Polynomial<T> PolynomialMultiply(const std::vector<Polynomial<T>>& v) {
  if (v.empty()) {
    return {};
  }
  auto ans = v[0];
  for (int i = 1; i < v.size(); i++) {
    ans = PolynomialMultiply(ans, v[i]);
  }
  return ans;
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

template<class T>
void Normalize(Polynomial<T>& a) {
  auto max = std::abs(*std::max_element(a.begin(), a.end(), [] (T a, T b) {
    return std::abs(a) < std::abs(b);
  })) / 1e6;
  for (auto& it: a) {
    it /= max;
  }
}

template<class T>
void SubtractPolynomial(Polynomial<T>& a, const Polynomial<T>& b) {
  if (a.size() < b.size()) {
    throw std::invalid_argument("");
  }
  int n = a.size();
  int m = b.size();
  for (int i = 0; i < b.size(); i++) {
    a[n - m + i] -= b[i];
  }
}

template<class T>
Polynomial<T> DividePolynomial(Polynomial<T> a,
                      const Polynomial<T>& b) {
  int n = a.size();
  int m = b.size();
  Polynomial<T> ans;
  for (int i = 0; n - i >= m; i++) {
    ans.push_back(a[i] / b[0]);
    Polynomial<T> x(n - m - i + 1);
    x[0] = ans.back();
    SubtractPolynomial(a, PolynomialMultiply(b, x));
  }
  return ans;
}
