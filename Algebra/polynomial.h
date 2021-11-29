#pragma once

#include <vector>

template<class T>
std::vector<T> PolynomialMultiply(const std::vector<T>& a,
                                  const std::vector<T>& b) {
  std::vector<T> c(a.size() + b.size() - 1);
  for (int i = 0; i < a.size(); i++) {
    for (int j = 0; j < b.size(); j++) {
      c[i + j] += a[i] * b[j];
    }
  }
  return c;
}

template<class T>
std::string PolynomialToString(const std::vector<T>& a,
                               const std::string& x = "x") {
  std::stringstream ss;
  for (int i = 0; i < a.size(); i++) {
    ss << a[i] << " * " << x << "^" << a.size() - i - 1 << " + ";
  }
  auto ans = ss.str();
  return ans.substr(0, ans.size() - 3) + "\n";
}
