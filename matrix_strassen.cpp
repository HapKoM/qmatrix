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
      data_[row_ctr] = new int8_t[(*it).size()];
      size_t col_ctr = 0;
      for (auto it1 = (*it).begin(); it1 != (*it).end(); ++it1) {
        data_[row_ctr][col_ctr] = *it1;
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
    data_[i] = new int8_t[col_];
  }
  clear();
}

Matrix::Matrix(const Matrix& other)
                     :row_(other.row_), col_(other.col_), data_(nullptr) {
  data_ = new int8_t*[row_];
  for (size_t i = 0; i < row_; ++i) {
    data_[i] = new int8_t[col_];
  }
  clear();
  for (size_t i = 0; i < row_; ++i) {
    memcpy(data_[i], other.data_[i], col_);
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
    memcpy(data_[i], rhs.data_[i], col_);
  }
  return *this;
}

bool Matrix::operator==(const Matrix& rhs) const {
  if (this == &rhs)
    return true;
  if ((row_ != rhs.row_) || (col_ != rhs.col_))
    return false;
  for (size_t i = 0; i < row_; ++i) {
     if (memcmp(data_[i], rhs.data_[i], col_) != 0)
      return false;
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
    for (size_t j = 0; j < col_; ++j) {
      m.data_[i][j] = (data_[i][j] + rhs.data_[i][j]) & mask;
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
    for (size_t j = 0; j < col_; ++j) {
      m.data_[i][j] = (data_[i][j] - rhs.data_[i][j]) & mask;
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
    new_data[i] = new int8_t[col_];
    memset(new_data[i], 0, col_);
    if (i < old_row) {
      if (old_col < col_) {
        // Increase old matrix, just copy old data to new
        memcpy(new_data[i], data_[i], old_col);
      } else {
        // Decrease old matrix, truncate old data to new size
        memcpy(new_data[i], data_[i], col_);
      }
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
    memset(data_[i], 0, col_);
  }
}

void Matrix::dump_size() const {
  std::cout << "[" << row_ << " x " << col_ << "]" << std::endl;
}

void Matrix::dump() const {
  for (size_t i = 0; i < row_; ++i) {
    for (size_t j = 0; j < col_; ++j) {
      std::cout << static_cast<int>(data_[i][j]) << " ";
    }
    std::cout << std::endl;
  }
}

char Matrix::get(size_t i, size_t j) const {
  if (i >= row_) {
    std::stringstream msg;
    msg << "Matrix::at(): Row number = " << i << "; should be less than " << row_;
    throw std::length_error(msg.str());
  }
  if (j >= col_) {
    std::stringstream msg;
    msg << "Matrix::at(): Column number = " << j << "; should be less than " << col_;
    throw std::length_error(msg.str());
  }
  return data_[i][j];
}

void Matrix::set(size_t i, size_t j, int8_t value) {
  if (i >= row_) {
    std::stringstream msg;
    msg << "Matrix::at(): Row number = " << i << "; should be less than " << row_;
    throw std::length_error(msg.str());
  }
  if (j >= col_) {
    std::stringstream msg;
    msg << "Matrix::at(): Column number = " << j << "; should be less than " << col_;
    throw std::length_error(msg.str());
  }
  data_[i][j] = value & mask;
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
      m.data_[j][i] = data_[i][j];
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
  for (size_t i = 0; i < lhs.row_; ++i) {
    for (size_t j = 0; j < rhs.col_; ++j) {
      int8_t sum = 0;
      for (size_t k = 0; k < lhs.col_; ++k) {
        sum += lhs.data_[i][k] * rhs.data_[k][j];
      }
      m.data_[i][j] = sum;
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
    memcpy(a_1_1.data_[i], a.data_[i], half_size);
    memcpy(a_1_2.data_[i], a.data_[i] + half_size, half_size);
    memcpy(a_2_1.data_[i], a.data_[i + half_size], half_size);
    memcpy(a_2_2.data_[i], a.data_[i + half_size] + half_size, half_size);
    memcpy(b_1_1.data_[i], b.data_[i], half_size);
    memcpy(b_1_2.data_[i], b.data_[i] + half_size, half_size);
    memcpy(b_2_1.data_[i], b.data_[i + half_size], half_size);
    memcpy(b_2_2.data_[i], b.data_[i + half_size] + half_size, half_size);
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
  for (size_t i = 0; i < half_size; ++i)
    memcpy(c.data_[i], c_1_1.data_[i], half_size);
  for (size_t i = 0; i < half_size; ++i)
    memcpy(c.data_[i] + half_size, c_1_2.data_[i], half_size);
  for (size_t i = 0; i < half_size; ++i)
    memcpy(c.data_[i + half_size], c_2_1.data_[i], half_size);
  for (size_t i = 0; i < half_size; ++i)
    memcpy(c.data_[i + half_size] + half_size, c_2_2.data_[i], half_size);
  c.resize(lhs.row_, rhs.col_);
  return c;
}
