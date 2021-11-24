#pragma once

#include <memory>

template<class T>
class Matrix {
 public:
  Matrix();
  Matrix(int n, int m);
  explicit Matrix(int n);

  std::pair<int, int> Size() const;

  int Rows() const;
  int Cols() const;

  Matrix<T> Copy();
  Matrix<T> TransposedCopy();

  static Matrix<T> Ones(int n, int m);
  static Matrix<T> Zeros(int n, int m);
  static Matrix<T> Random(int min, int max, int seed = time(nullptr));

  std::string ToWolframString() const;

  T& operator()(int i);  // For vectors only
  const T& operator()(int i) const;  // For vectors only
  T& operator()(int i, int j);
  const T& operator()(int i, int j) const;

  Matrix<T> SumMatrix(int i, int j, int n, int m) const;

  Matrix<T> Row(int i) const;
  Matrix<T> Col(int j) const;

  friend Matrix<T> operator+(const Matrix<T>& a, const Matrix<T>& b);
  friend Matrix<T> operator-(const Matrix<T>& a, const Matrix<T>& b);
  friend Matrix<T> operator*(const Matrix<T>& a, const Matrix<T>& b);
  friend Matrix<T> operator/(const Matrix<T>& a, const Matrix<T>& b);
  friend T ScalarProduct(const Matrix<T>& a, const Matrix<T>& b);
  friend Matrix<T> operator*(const Matrix<T>& a, T b);
  friend Matrix<T> operator/(const Matrix<T>& a, T b);

  Matrix<T>& operator+=(const Matrix<T>& a);
  Matrix<T>& operator-=(const Matrix<T>& a);
  Matrix<T>& operator*=(const Matrix<T>& a);
  Matrix<T>& operator/=(const Matrix<T>& a);
  Matrix<T> operator*=(T b);
  Matrix<T> operator/=(T b);

  friend operator==(const Matrix<T>& a, const Matrix<T>& b);
  friend operator!=(const Matrix<T>& a, const Matrix<T>& b);

 private:
  std::shared_ptr<T[]> data_;
  int data_rows_;
  int data_cols_;
  int cols_;
  int rows_;

  bool IsSubMatrix() const;
};

template<class U>
std::ostream& operator<<(std::ostream& stream,
                         const Matrix<U>& matrix) {
  size_t maxlen = 0;
  for (size_t i = 0; i < matrix.Rows(); i++) {
    for (size_t j = 0; j < matrix.Cols(); j++) {
      std::stringstream ss;
      ss << std::fixed << std::setprecision(5) << matrix.At(i, j);
      maxlen = std::max(maxlen, ss.str().length());
    }
  }
  stream << "[";
  for (size_t i = 0; i < matrix.Rows(); i++) {
    if (i != 0) {
      stream << ' ';
    }
    for (size_t j = 0; j < matrix.Cols(); j++) {
      stream << std::fixed << std::setprecision(5) << std::setw(maxlen)
             << matrix.At(i, j);
      if (i + 1 < matrix.Rows() || j + 1 < matrix.Cols()) {
        stream << ", ";
      }
    }
    if (i + 1 < matrix.Rows()) {
      stream << '\n';
    }
  }
  stream << "]\n";
  return stream;
}
