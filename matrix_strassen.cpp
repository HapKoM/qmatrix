#include <stdexcept>
#include <sstream>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <future>
#include <functional>

#include "matrix_strassen.h"

#ifndef STRASSEN_MATRIX_SIZE
#define STRASSEN_MATRIX_SIZE 64
#endif

inline size_t packed_bytes_size(size_t col) {
  size_t n_bytes = col / 4;
  if (col % 4)
    n_bytes++;
  return n_bytes;
}

#ifndef TEST_MODE
inline int8_t packed_sum(int8_t a, int8_t b)
#else
int8_t Matrix::packed_sum(int8_t a, int8_t b)
#endif
{
  int8_t v0 = (a + b) & 0x03;
  int8_t v1 = (((a >> 2) + (b >> 2)) & 0x03) << 2;
  int8_t v2 = (((a >> 4) + (b >> 4)) & 0x03) << 4;
  int8_t v3 = (((a >> 6) + (b >> 6)) & 0x03) << 6;
  return v3 | v2 | v1 | v0;
}

#ifndef TEST_MODE
inline int8_t packed_diff(int8_t a, int8_t b)
#else
int8_t Matrix::packed_diff(int8_t a, int8_t b)
#endif
{
  int8_t v0 = (a - b) & 0x03;
  int8_t v1 = (((a >> 2) - (b >> 2)) & 0x03) << 2;
  int8_t v2 = (((a >> 4) - (b >> 4)) & 0x03) << 4;
  int8_t v3 = (((a >> 6) - (b >> 6)) & 0x03) << 6;
  return v3 | v2 | v1 | v0;
}

#ifndef TEST_MODE
inline int8_t packed_multiply(int8_t a, int8_t b)
#else
int8_t Matrix::packed_multiply(int8_t a, int8_t b)
#endif
{
  int8_t v0 = (a*b) & 0x03;
  int8_t v1 = (((a >> 2) * (b >> 2)) & 0x03) << 2;
  int8_t v2 = (((a >> 4) * (b >> 4)) & 0x03) << 4;
  int8_t v3 = (((a >> 6) * (b >> 6)) & 0x03) << 6;
  return v3 | v2 | v1 | v0;
}

Matrix::Matrix(std::initializer_list<std::initializer_list<int8_t> > data)
       :row_(data.size()), col_(0), data_(nullptr) {
  if (row_ > 0) {
    col_ = (*data.begin()).size();
    size_t row_ctr = 0;
    // First pass to check equality of rows sizes
    for (auto it = data.begin(); it != data.end(); ++it) {
      if (col_ != (*it).size()) {
        std::stringstream msg;
        msg << "Matrix::Matrix(): Different column numbers in initializer_list rows";
        throw std::length_error(msg.str());
      }
    }
    data_ = new int8_t*[row_];
    for (auto it = data.begin(); it != data.end(); ++it) {
      size_t packed_size = packed_bytes_size((*it).size());
      data_[row_ctr] = new int8_t[packed_size];
      memset(data_[row_ctr], 0, packed_size);
      size_t col_ctr = 0;
      for (auto it1 = (*it).begin(); it1 != (*it).end(); ++it1) {
        set(row_ctr, col_ctr, *it1);
        col_ctr++;
      }
      row_ctr++;
    }
  }
}

Matrix::Matrix(size_t row, size_t col)
                     :row_(row), col_(col), data_(nullptr) {
  data_ = new int8_t*[row_];
  for (size_t i = 0; i < row_; ++i) {
    size_t packed_size = packed_bytes_size(col_);
    data_[i] = new int8_t[packed_size];
    memset(data_[i], 0, packed_size);
  }
  clear();
}

Matrix::Matrix(const Matrix& other)
                     :row_(other.row_), col_(other.col_), data_(nullptr) {
  data_ = new int8_t*[row_];
  for (size_t i = 0; i < row_; ++i) {
    size_t packed_size = packed_bytes_size(col_);
    data_[i] = new int8_t[packed_size];
    memset(data_[i], 0, packed_size);
  }
  clear();
  for (size_t i = 0; i < row_; ++i) {
    memcpy(data_[i], other.data_[i], packed_bytes_size(col_));
  }
}

Matrix::~Matrix() {
  for (size_t i = 0; i < row_; ++i) {
    if (data_[i])
      delete[] data_[i];
  }
  if (data_)
    delete[] data_;
}


Matrix& Matrix::operator=(const Matrix& rhs) {
  if (this == &rhs)
    return *this;
  clear();
  row_ = rhs.row_;
  col_ = rhs.col_;
  for (size_t i = 0; i < row_; ++i) {
    memcpy(data_[i], rhs.data_[i], packed_bytes_size(col_));
  }
  return *this;
}

bool Matrix::operator==(const Matrix& rhs) const {
  if (this == &rhs)
    return true;
  if ((row_ != rhs.row_) || (col_ != rhs.col_))
    return false;
  for (size_t i = 0; i < row_; ++i) {
    size_t packed_size = packed_bytes_size(col_);
    if (memcmp(data_[i], rhs.data_[i], packed_size) != 0) {
      return false;
    }
  }
  return true;
}

Matrix Matrix::operator*(const Matrix& rhs) const {
  if (col_ != rhs.row_) {
    std::stringstream msg;
    msg << "Matrix::operator*: Column number of first matrix should be equal to row number of the second matrix ("
        << col_ << " and " << rhs.row_ << " provided)";
    throw std::length_error(msg.str());
  }
  size_t max_size = std::max(std::max(col_, row_), std::max(rhs.col_, rhs.row_));
#ifdef TRIVIAL_ALGORITHM
  return multiply_trivial(*this, rhs);
#else
  /*
   * Strassen algorithm is effective for matrices of size
   * more than 64x64. For small matrices we use the trivial algorithm.
   * Note: If STRASSEN_MATRIX_SIZE is too low and matrices is big (256*256 and more)
   * multiply_strassen may throw std::system_error
   */
  if (max_size <= STRASSEN_MATRIX_SIZE) {
    return multiply_trivial(*this, rhs);
  } else {
    return multiply_strassen(*this, rhs);
  }
#endif
}

Matrix Matrix::operator+(const Matrix& rhs) const {
  if ((row_ != rhs.row_) || (col_ != rhs.col_)) {
    std::stringstream msg;
    msg << "Matrix::operator+: Sizes of matricies should be equal("
        << row_ << "x" << col_ << " and " << rhs.row_ << "x" << rhs.col_ << " provided)";
    throw std::length_error(msg.str());
  }
  Matrix m(row_, col_);
  for (size_t i = 0; i < row_; ++i) {
    size_t n_bytes = packed_bytes_size(col_);
    for (size_t j = 0; j < n_bytes; ++j) {
      m.data_[i][j] = packed_sum(data_[i][j], rhs.data_[i][j]);
    }
  }
  return m;
}

Matrix Matrix::operator-(const Matrix& rhs) const {
  if ((this->row() != rhs.row_) || (this->col() != rhs.col_)) {
    std::stringstream msg;
    msg << "Matrix::operator-: Sizes of matricies should be equal("
        << this->col() << "x" << this->row() << " and " << rhs.col_ << "x" << rhs.row_ << " provided)";
    throw std::length_error(msg.str());
  }
  Matrix m(this->row(), this->col());
  for (size_t i = 0; i < row_; ++i) {
    size_t n_bytes = packed_bytes_size(col_);
    for (size_t j = 0; j < n_bytes; ++j) {
      m.data_[i][j] = packed_diff(data_[i][j], rhs.data_[i][j]);
    }
  }
  return m;
}

void Matrix::resize(size_t row, size_t col) {
  size_t old_col = col_;
  size_t old_row = row_;
  col_ = col;
  row_ = row;
  int8_t** new_data = new int8_t*[row_];
  for (size_t i = 0; i < col_; ++i) {
    new_data[i] = new int8_t[packed_bytes_size(col_)];
    memset(new_data[i], 0, packed_bytes_size(col_));
    if (i < old_row) {
      size_t sz = std::min(old_col, col_);
      size_t packed_sz = packed_bytes_size(sz);
      memcpy(new_data[i], data_[i], packed_sz);
      // If column size is not multiple of 4, we need to zeroize
      // 2, 4 or 6 most significant bits in the last byte of current row
      // Number of bits to be zeroized depends on column size modulo 4
      int8_t shift = sz % 4;
      int8_t mask = 0xFF >> shift*2;
      new_data[i][packed_sz - 1] &= mask;
    }
  }
  for (size_t i = 0; i < old_row; ++i) {
    delete[] data_[i];
  }
  delete[] data_;
  data_ = new_data;
}

void Matrix::clear() {
  for (size_t i = 0; i < row_; ++i) {
    memset(data_[i], 0, packed_bytes_size(col_));
  }
}

void Matrix::dump_size() const {
  std::cout << "[" << row_ << " x " << col_ << "]" << std::endl;
}

void Matrix::dump_raw_bytes() const {
  for (size_t i = 0; i < row_; ++i) {
    size_t n_bytes = packed_bytes_size(col_);
    for (size_t j = 0; j < n_bytes; ++j) {
      printf("%02X ", data_[i][j]);
    }
    std::cout << std::endl;
  }
}

void Matrix::dump() const {
  std::cout << "[";
  for (size_t i = 0; i < row_; ++i) {
    for (size_t j = 0; j < col_; ++j) {
      std::cout << static_cast<int>(get(i, j)) << " ";
    }
    if (i != row_ - 1)
      std::cout << "; ";
  }
  std::cout << "];" << std::endl;
}

size_t Matrix::row() const {
  return row_;
}

size_t Matrix::col() const {
  return col_;
}

Matrix Matrix::transposed() const {
  Matrix m(col_, row_);
  for (size_t i = 0; i < row_; ++i) {
    for (size_t j = 0; j < col_; ++j) {
      m.set(j, i, get(i, j));
    }
  }
  return m;
}

Matrix Matrix::multiply_trivial(const Matrix& lhs, const Matrix& rhs) {
  if (lhs.col_ != rhs.row_) {
    std::stringstream msg;
    msg << "Matrix::multiply_trivial: Column number of first marrix should be equal to row number of the second matrix ("
        << lhs.col_ << " and " << rhs.row_ << " provided)";
    throw std::length_error(msg.str());
  }

  Matrix m(lhs.row_, rhs.col_);
  Matrix rhs_tr(rhs.transposed());
  size_t l_bytes = packed_bytes_size(lhs.col_);
  for (size_t i = 0; i < lhs.row_; ++i) {
    for (size_t j = 0; j < rhs_tr.row_; ++j) {
      int8_t sum = 0;
      for (size_t k = 0; k < l_bytes; ++k) {
        int8_t prod = packed_multiply(lhs.data_[i][k], rhs_tr.data_[j][k]);
        sum = packed_sum(sum, prod);
      }
      int8_t sum_pack = ((sum & 0x03) + ((sum >> 2) & 0x03) + ((sum >> 4) & 0x03) + ((sum >> 6) & 0x03)) & 0x03;
      m.set(i, j, sum_pack);
    }
  }
  return m;
}

Matrix Matrix::calculate_p1(const Matrix& a11, const Matrix& a22, const Matrix& b11, const Matrix& b22) {
  return (a11 + a22)*(b11 + b22);
}

Matrix Matrix::calculate_p2(const Matrix& a21, const Matrix& a22, const Matrix& b11) {
  return (a21 + a22) * b11;
}

Matrix Matrix::calculate_p3(const Matrix& a11, const Matrix& b12, const Matrix& b22) {
  return a11 * (b12 - b22);
}

Matrix Matrix::calculate_p4(const Matrix& a22, const Matrix& b21, const Matrix& b11) {
  return a22 * (b21 - b11);
}

Matrix Matrix::calculate_p5(const Matrix& a11, const Matrix& a12, const Matrix& b22) {
  return (a11 + a12) * b22;
}

Matrix Matrix::calculate_p6(const Matrix& a21, const Matrix& a11, const Matrix& b11, const Matrix& b12) {
  return (a21 - a11) * (b11 + b12);
}

Matrix Matrix::calculate_p7(const Matrix& a12, const Matrix& a22, const Matrix& b21, const Matrix& b22) {
  return (a12 - a22) * (b21 + b22);
}

Matrix Matrix::multiply_strassen(const Matrix& lhs, const Matrix& rhs) {
  /*
   * Strassen algorithm implementation
   * See https://en.wikipedia.org/wiki/Strassen_algorithm for details
   */
  if ((lhs.row_ == 1) && (lhs.col_ == 1) && (rhs.col_ == 1) && (rhs.col_ == 1)) {
    Matrix m(1, 1);
    m.data_[0][0] = lhs.data_[0][0] * rhs.data_[0][0];
    return m;
  }
  size_t max_size = std::max(std::max(lhs.col_, lhs.row_), std::max(rhs.col_, rhs.row_));
  size_t power = 1;
  // Find next highest power of 2
  while (max_size > power) power *= 2;
  // Prepare matrices to be squared with sizes [power x power]
  Matrix a(lhs);
  Matrix b(rhs);
  a.resize(power, power);
  b.resize(power, power);
  size_t half_size = power / 2;
  //size_t packed_size = packed_bytes_size(half_size);
  // Initialize and fill quarters of matrices
  Matrix a_1_1(half_size, half_size);
  Matrix a_1_2(half_size, half_size);
  Matrix a_2_1(half_size, half_size);
  Matrix a_2_2(half_size, half_size);
  Matrix b_1_1(half_size, half_size);
  Matrix b_1_2(half_size, half_size);
  Matrix b_2_1(half_size, half_size);
  Matrix b_2_2(half_size, half_size);
  for (size_t i = 0; i < half_size; ++i) {
    for (size_t j = 0; j < half_size; ++j) {
      a_1_1.set(i, j, a.get(i, j));
      a_1_2.set(i, j, a.get(i, j + half_size));
      a_2_1.set(i, j, a.get(i + half_size, j));
      a_2_2.set(i, j, a.get(i + half_size, j + half_size));
      b_1_1.set(i, j, b.get(i, j));
      b_1_2.set(i, j, b.get(i, j + half_size));
      b_2_1.set(i, j, b.get(i + half_size, j));
      b_2_2.set(i, j, b.get(i + half_size, j + half_size));
    }
  }
#ifdef PARALLEL_STRASSEN
#pragma message "parellelized Strassen algorithm implementation"
#pragma message "Undefine PARALLEL_STRASSEN to disable parallelization"
  std::future<Matrix> f1 = std::async(std::launch::async, Matrix::calculate_p1, std::cref(a_1_1), std::cref(a_2_2), std::cref(b_1_1), std::cref(b_2_2));
  std::future<Matrix> f2 = std::async(std::launch::async, Matrix::calculate_p2, std::cref(a_2_1), std::cref(a_2_2), std::cref(b_1_1));
  std::future<Matrix> f3 = std::async(std::launch::async, Matrix::calculate_p3, std::cref(a_1_1), std::cref(b_1_2), std::cref(b_2_2));
  std::future<Matrix> f4 = std::async(std::launch::async, Matrix::calculate_p4, std::cref(a_2_2), std::cref(b_2_1), std::cref(b_1_1));
  std::future<Matrix> f5 = std::async(std::launch::async, Matrix::calculate_p5, std::cref(a_1_1), std::cref(a_1_2), std::cref(b_2_2));
  std::future<Matrix> f6 = std::async(std::launch::async, Matrix::calculate_p6, std::cref(a_2_1), std::cref(a_1_1), std::cref(b_1_1), std::cref(b_1_2));
  std::future<Matrix> f7 = std::async(std::launch::async, Matrix::calculate_p7, std::cref(a_1_2), std::cref(a_2_2), std::cref(b_2_1), std::cref(b_2_2));
  f1.wait();
  f2.wait();
  f3.wait();
  f4.wait();
  f5.wait();
  f6.wait();
  f7.wait();
  Matrix p_1(f1.get());
  Matrix p_2(f2.get());
  Matrix p_3(f3.get());
  Matrix p_4(f4.get());
  Matrix p_5(f5.get());
  Matrix p_6(f6.get());
  Matrix p_7(f7.get());
#else
#pragma message "Non-parellelized Strassen algorithm implementation"
#pragma message "Define PARALLEL_STRASSEN to enable parallelization"
  Matrix p_1((a_1_1 + a_2_2) * (b_1_1 + b_2_2));
  Matrix p_2((a_2_1 + a_2_2) * b_1_1);
  Matrix p_3(a_1_1 * (b_1_2 - b_2_2));
  Matrix p_4(a_2_2 * (b_2_1 - b_1_1));
  Matrix p_5((a_1_1 + a_1_2) * b_2_2);
  Matrix p_6((a_2_1 - a_1_1) * (b_1_1 + b_1_2));
  Matrix p_7((a_1_2 - a_2_2) * (b_2_1 + b_2_2));
#endif

  Matrix c_1_1(p_1 + p_4 - p_5 + p_7);
  Matrix c_1_2(p_3 + p_5);
  Matrix c_2_1(p_2 + p_4);
  Matrix c_2_2(p_1 - p_2 + p_3 + p_6);

  Matrix c(power, power);

  for (size_t i = 0; i < half_size; ++i) {
    for (size_t j = 0; j < half_size; ++j) {
      c.set(i, j, c_1_1.get(i, j));
      c.set(i, j + half_size, c_1_2.get(i, j));
      c.set(i + half_size, j, c_2_1.get(i, j));
      c.set(i + half_size, j + half_size, c_2_2.get(i, j));
    }
  }
  c.resize(lhs.row_, rhs.col_);
  return c;
}
