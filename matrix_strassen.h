#ifndef MATRIX_STRASSEN_H
#define MATRIX_STRASSEN_H

#include <cstddef>

#include <stdint.h>

class Matrix {
public:
  /**
   * Creates Matrix object and fills it with provided data:
   * Matrix m({ {1, 2, 3}, {4, 5, 6} });
   * Throws std::length_error if some of rows in data have different size
   */
  Matrix(std::initializer_list<std::initializer_list<int8_t> > data);
  Matrix(size_t row, size_t col);
  Matrix(const Matrix& other);
  ~Matrix();
  Matrix& operator=(const Matrix& rhs);
  bool operator==(const Matrix& rhs) const;
  /**
   * Matrix multiplication
   * Throws std::length_error if column number of this matrix
   * is not equal to row number of the rhs
   */
  Matrix operator*(const Matrix& rhs) const;
  /**
   * Matrix addition
   * Sizes of matrices should be equal, otherwise std::length_error is thrown
   */
  Matrix operator+(const Matrix& rhs) const;
  /**
   * Matrix subtraction
   * Sizes of matrices should be equal, otherwise std::length_error is thrown
   */
  Matrix operator-(const Matrix& rhs) const;
  /**
   * Set matrix size to [row, col]
   * If new size is less than old, matrix data will be truncated to fit new size;
   * Otherwise matrix data will be increased and new cells will be filled with zeroes
   */
  void resize(size_t row, size_t col);
  //! Zeroize data
  void clear();
  void dump_size() const;
  void dump_raw_bytes() const;
  void dump() const;
  /**
   * Get element in i-th row and j-th column
   * i should be [0..row), j should be [0..col)
   */
  inline int8_t get(size_t i, size_t j) const {
    size_t n_byte = j / 4;
    size_t shift = (j % 4) * 2;
    return (data_[i][n_byte] >> shift) & 0x03;
  }

  /**
   * Set element in i-th row and j-th column
   * i should be [0..row), j should be [0..col)
   */
  inline void set(size_t i, size_t j, int8_t value) {
    size_t n_byte = j / 4;
    size_t shift = (j % 4) * 2;
    int8_t old_byte = data_[i][n_byte];
    int8_t mask = ~(0x03 << shift);
    value = (value & 0x03) << shift;
    data_[i][n_byte] = (old_byte & mask) | value;
  }

  size_t row() const;
  size_t col() const;
  Matrix transposed() const;
  static Matrix multiply_trivial(const Matrix& lhs, const Matrix& rhs);
  static Matrix multiply_strassen(const Matrix& lhs, const Matrix& rhs);
#ifdef TEST_MODE
  /*
   * For the tesing this functions declared as static methods.
   * When testing is disabled this methods defined as inline functions
   * and visible only inside matrix_strassen.cpp file
   */
  static int8_t packed_sum(int8_t a, int8_t b);
  static int8_t packed_diff(int8_t a, int8_t b);
  static int8_t packed_multiply(int8_t a, int8_t b);
#endif
private:
  static Matrix calculate_p1(const Matrix& a11, const Matrix& a22, const Matrix& b11, const Matrix& b22);
  static Matrix calculate_p2(const Matrix& a21, const Matrix& a22, const Matrix& b11);
  static Matrix calculate_p3(const Matrix& a11, const Matrix& b12, const Matrix& b22);
  static Matrix calculate_p4(const Matrix& a22, const Matrix& b21, const Matrix& b11);
  static Matrix calculate_p5(const Matrix& a11, const Matrix& a12, const Matrix& b22);
  static Matrix calculate_p6(const Matrix& a21, const Matrix& a11, const Matrix& b11, const Matrix& b12);
  static Matrix calculate_p7(const Matrix& a12, const Matrix& a22, const Matrix& b21, const Matrix& b22);
  //! Rows number
  size_t row_;
  //! Columns number
  size_t col_;
  //! Pointer to matrix data
  int8_t** data_;
};

#endif // MATRIX_STRASSEN_H
