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

  Matrix<T> Transposed();

  static Matrix<T> Ones(int n, int m);
  static Matrix<T> Zeros(int n, int m);
  Matrix<T> Matrix<T>::Random(int n,
                              int m,
                              T min,
                              T max,
                              int seed = time(nullptr));

  std::string ToWolframString() const;

  T& operator()(int i);  // For vectors only
  const T& operator()(int i) const;  // For vectors only
  T& operator()(int i, int j);
  const T& operator()(int i, int j) const;

  T& At(int i);  // For vectors only
  const T& At(int i) const;  // For vectors only
  T& At(int i, int j);
  const T& At(int i, int j) const;

  const Matrix<T> SubMatrix(int i, int j, int n, int m) const;
  Matrix<T> SubMatrix(int i, int j, int n, int m);

  bool IsRowVector() const;
  bool IsColVector() const;

  const Matrix<T> Row(int i) const;
  Matrix<T> Row(int i);
  const Matrix<T> Col(int j) const;
  Matrix<T> Col(int j);

  friend Matrix<T>& operator=(Matrix<T>& a, const Matrix<T>& b);
  friend Matrix<T>& operator=(Matrix<T>& a, Matrix<T>&& b);

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

template<class T>
Matrix<T>::Matrix() : Matrix<T>(0, 0) {}

template<class T>
Matrix<T>::Matrix(int n, int m) : data_(new T[n * m]) {
  data_rows_ = rows_ = n;
  data_cols_ = cols_ = m;
  for (int i = 0; i < this->Rows(); i++) {
    for (int j = 0; j < this->Cols(); j++) {
      (*this)(i, j) = T();
    }
  }
}

template<class T>
Matrix<T>::Matrix(int n) : Matrix<T>(n, n) {}

template<class T>
std::pair<int, int> Matrix<T>::Size() const {
  return {rows_, cols_};
}

template<class T>
int Matrix<T>::Rows() const {
  return rows_;
}

template<class T>
int Matrix<T>::Cols() const {
  return cols_;
}

template<class T>
Matrix<T> Matrix<T>::Transposed() {
  return Matrix<T>();
}

template<class T>
Matrix<T> Matrix<T>::Ones(int n, int m) {
  Matrix<T> a(n, m);
  for (int i = 0; i < std::min(n, m); i++) {
    a(i, i) = 1;
  }
  return a;
}

template<class T>
Matrix<T> Matrix<T>::Zeros(int n, int m) {
  return Matrix<T>(n, m);
}

template<class T>
Matrix<T> Matrix<T>::Random(int n, int m, T min, T max, int seed) {
  Matrix<T> a(n, m);
  std::mt19937 gen(seed);
  std::uniform_real_distribution<T> dist(min, max);
  for (int i = 0; i < this->Rows(); i++) {
    for (int j = 0; j < this->Cols(); j++) {
      a(i, j) = dist(gen);
    }
  }
  return a;
}

template<class T>
std::string Matrix<T>::ToWolframString() const {
  std::stringstream res;
  res << "{";
  for (int i = 0; i < this->Rows(); i++) {
    res << "{";
    for (int j = 0; j < this->Cols(); j++) {
      res << std::fixed << std::setprecision(4) << this->At(i, j);
      if (j + 1 != this->Cols()) {
        res << ",";
      }
    }
    res << "}";
    if (i + 1 != this->Rows()) {
      res << ",";
    }
  }
  res << "}";
  return res.str();
}

template<class T>
T& Matrix<T>::operator()(int i) {
  return this->At(i);
}

template<class T>
const T& Matrix<T>::operator()(int i) const {
  return this->At(i);
}

template<class T>
T& Matrix<T>::operator()(int i, int j) {
  return this->At(i, j);
}

template<class T>
const T& Matrix<T>::operator()(int i, int j) const {
  return this->At(i, j);
}

template<class T>
const Matrix<T> Matrix<T>::SubMatrix(int i, int j, int n, int m) const {
  return Matrix<T>();
}
template<class T>
Matrix<T> Matrix<T>::SubMatrix(int i, int j, int n, int m) {
  return Matrix<T>();
}
template<class T>
const Matrix<T> Matrix<T>::Row(int i) const {
  return Matrix<T>();
}
template<class T>
Matrix<T> Matrix<T>::Row(int i) {
  return Matrix<T>();
}
template<class T>
const Matrix<T> Matrix<T>::Col(int j) const {
  return Matrix<T>();
}
template<class T>
Matrix<T> Matrix<T>::Col(int j) {
  return Matrix<T>();
}
template<class T>
Matrix<T>& Matrix<T>::operator+=(const Matrix<T>& a) {
  return <#initializer#>;
}
template<class T>
Matrix<T>& Matrix<T>::operator-=(const Matrix<T>& a) {
  return <#initializer#>;
}
template<class T>
Matrix<T>& Matrix<T>::operator*=(const Matrix<T>& a) {
  return <#initializer#>;
}
template<class T>
Matrix<T>& Matrix<T>::operator/=(const Matrix<T>& a) {
  return <#initializer#>;
}
template<class T>
Matrix<T> Matrix<T>::operator*=(T b) {
  return Matrix<T>();
}
template<class T>
Matrix<T> Matrix<T>::operator/=(T b) {
  return Matrix<T>();
}
template<class T>
bool Matrix<T>::IsSubMatrix() const {
  return false;
}

template<class T>
bool Matrix<T>::IsRowVector() const {
  return false;
}
template<class T>
bool Matrix<T>::IsColVector() const {
  return false;
}
template<class T>
T& Matrix<T>::At(int i) {
  return <#initializer#>;
}
template<class T>
const T& Matrix<T>::At(int i) const {
  return <#initializer#>;
}
template<class T>
T& Matrix<T>::At(int i, int j) {
  return <#initializer#>;
}
template<class T>
const T& Matrix<T>::At(int i, int j) const {
  return <#initializer#>;
}

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
