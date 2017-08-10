#include <stdexcept>
#include <random>

#include <gtest/gtest.h>

#include "matrix_strassen.h"

class MatrixTest : public ::testing::Test
{
protected:
  void SetUp()
  {
    srand (time(NULL));
  }

  void TearDown()
  {
  }
};

TEST(MatrixTest, MatrixEqualityTest) {
  Matrix a({{1, 2, 3}, {4, 5, 6}});
  Matrix b(2, 3);
  b.set(0, 0, 1);
  b.set(0, 1, 2);
  b.set(0, 2, 3);
  b.set(1, 0, 4);
  b.set(1, 1, 5);
  b.set(1, 2, 6);
  Matrix c = a;
  Matrix d(c);
  ASSERT_EQ(a, b);
  ASSERT_EQ(a, c);
  ASSERT_EQ(c, d);
}

TEST(MatrixTest, MatrixTranspositionTest) {
  Matrix a({{1, 2, 3}, {4, 5, 6}});
  Matrix b({{1, 4}, {2, 5}, {3, 6}});
  ASSERT_EQ(a.transposed(), b);
}

TEST(MatrixTest, MatrixClearTest) {
  Matrix a({{1, 2, 3}, {4, 5, 6}});
  Matrix b(2, 3);
  Matrix c({{0, 0, 0}, {0, 0, 0}});
  a.clear();
  ASSERT_EQ(a, b);
  ASSERT_EQ(a, c);
}

TEST(MatrixTest, MatrixResizeTest) {
  Matrix a({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});
  a.resize(2, 2);
  Matrix b({{1, 2}, {4, 5}});
  ASSERT_EQ(a, b);
  b.resize(3, 3);
  Matrix c({{1, 2, 0}, {4, 5, 0}, {0, 0, 0}});
  ASSERT_EQ(b, c);
}

TEST(MatrixTest, MatrixErrorTest) {
  Matrix a({{1, 2, 3}, {4, 5, 6}});
  Matrix b = a.transposed();
  ASSERT_THROW(a*a, std::length_error);
  ASSERT_THROW(a + b, std::length_error);
  ASSERT_THROW(a - b, std::length_error);
  ASSERT_THROW(Matrix({{1, 2, 3}, {1, 2}}), std::length_error);
}

TEST(MatrixTest, PackedSumTest) {
  int8_t a0 = 0x00;
  int8_t a1 = 0x01;
  int8_t a2 = 0x02;
  int8_t a3 = 0x03;
  for (uint8_t i = 0; i < 0xFF; ++i) {
    for (uint8_t j = 0; j < 0xFF; ++j) {
      ASSERT_EQ(Matrix::packed_sum(i, j), Matrix::packed_sum(j, i));
    }
  }
  ASSERT_EQ(a0, Matrix::packed_sum(a0, a0));
  ASSERT_EQ(a1, Matrix::packed_sum(a0, a1));
  ASSERT_EQ(a2, Matrix::packed_sum(a0, a2));
  ASSERT_EQ(a3, Matrix::packed_sum(a0, a3));

  ASSERT_EQ(a1, Matrix::packed_sum(a1, a0));
  ASSERT_EQ(a2, Matrix::packed_sum(a1, a1));
  ASSERT_EQ(a3, Matrix::packed_sum(a1, a2));
  ASSERT_EQ(a0, Matrix::packed_sum(a1, a3));

  ASSERT_EQ(a2, Matrix::packed_sum(a2, a0));
  ASSERT_EQ(a3, Matrix::packed_sum(a2, a1));
  ASSERT_EQ(a0, Matrix::packed_sum(a2, a2));
  ASSERT_EQ(a1, Matrix::packed_sum(a2, a3));

  ASSERT_EQ(a3, Matrix::packed_sum(a3, a0));
  ASSERT_EQ(a0, Matrix::packed_sum(a3, a1));
  ASSERT_EQ(a1, Matrix::packed_sum(a3, a2));
  ASSERT_EQ(a2, Matrix::packed_sum(a3, a3));

  int8_t b0 = a0 | a0 << 2 | a0 << 4 | a0 << 6;
  int8_t b1 = a1 | a1 << 2 | a1 << 4 | a1 << 6;
  int8_t b2 = a2 | a2 << 2 | a2 << 4 | a2 << 6;
  int8_t b3 = a3 | a3 << 2 | a3 << 4 | a3 << 6;

  ASSERT_EQ(b0, Matrix::packed_sum(b0, b0));
  ASSERT_EQ(b1, Matrix::packed_sum(b0, b1));
  ASSERT_EQ(b2, Matrix::packed_sum(b0, b2));
  ASSERT_EQ(b3, Matrix::packed_sum(b0, b3));

  ASSERT_EQ(b1, Matrix::packed_sum(b1, b0));
  ASSERT_EQ(b2, Matrix::packed_sum(b1, b1));
  ASSERT_EQ(b3, Matrix::packed_sum(b1, b2));
  ASSERT_EQ(b0, Matrix::packed_sum(b1, b3));

  ASSERT_EQ(b2, Matrix::packed_sum(b2, b0));
  ASSERT_EQ(b3, Matrix::packed_sum(b2, b1));
  ASSERT_EQ(b0, Matrix::packed_sum(b2, b2));
  ASSERT_EQ(b1, Matrix::packed_sum(b2, b3));

  ASSERT_EQ(b3, Matrix::packed_sum(b3, b0));
  ASSERT_EQ(b0, Matrix::packed_sum(b3, b1));
  ASSERT_EQ(b1, Matrix::packed_sum(b3, b2));
  ASSERT_EQ(b2, Matrix::packed_sum(b3, b3));
}

TEST(MatrixTest, PackedDiffTest) {
  int8_t a0 = 0x00;
  int8_t a1 = 0x01;
  int8_t a2 = 0x02;
  int8_t a3 = 0x03;

  ASSERT_EQ(a0, Matrix::packed_diff(a0, a0));
  ASSERT_EQ(a3, Matrix::packed_diff(a0, a1));
  ASSERT_EQ(a2, Matrix::packed_diff(a0, a2));
  ASSERT_EQ(a1, Matrix::packed_diff(a0, a3));

  ASSERT_EQ(a1, Matrix::packed_diff(a1, a0));
  ASSERT_EQ(a0, Matrix::packed_diff(a1, a1));
  ASSERT_EQ(a3, Matrix::packed_diff(a1, a2));
  ASSERT_EQ(a2, Matrix::packed_diff(a1, a3));

  ASSERT_EQ(a2, Matrix::packed_diff(a2, a0));
  ASSERT_EQ(a1, Matrix::packed_diff(a2, a1));
  ASSERT_EQ(a0, Matrix::packed_diff(a2, a2));
  ASSERT_EQ(a3, Matrix::packed_diff(a2, a3));

  ASSERT_EQ(a3, Matrix::packed_diff(a3, a0));
  ASSERT_EQ(a2, Matrix::packed_diff(a3, a1));
  ASSERT_EQ(a1, Matrix::packed_diff(a3, a2));
  ASSERT_EQ(a0, Matrix::packed_diff(a3, a3));

  int8_t b0 = a0 | a0 << 2 | a0 << 4 | a0 << 6;
  int8_t b1 = a1 | a1 << 2 | a1 << 4 | a1 << 6;
  int8_t b2 = a2 | a2 << 2 | a2 << 4 | a2 << 6;
  int8_t b3 = a3 | a3 << 2 | a3 << 4 | a3 << 6;

  ASSERT_EQ(b0, Matrix::packed_diff(b0, b0));
  ASSERT_EQ(b3, Matrix::packed_diff(b0, b1));
  ASSERT_EQ(b2, Matrix::packed_diff(b0, b2));
  ASSERT_EQ(b1, Matrix::packed_diff(b0, b3));

  ASSERT_EQ(b1, Matrix::packed_diff(b1, b0));
  ASSERT_EQ(b0, Matrix::packed_diff(b1, b1));
  ASSERT_EQ(b3, Matrix::packed_diff(b1, b2));
  ASSERT_EQ(b2, Matrix::packed_diff(b1, b3));

  ASSERT_EQ(b2, Matrix::packed_diff(b2, b0));
  ASSERT_EQ(b1, Matrix::packed_diff(b2, b1));
  ASSERT_EQ(b0, Matrix::packed_diff(b2, b2));
  ASSERT_EQ(b3, Matrix::packed_diff(b2, b3));

  ASSERT_EQ(b3, Matrix::packed_diff(b3, b0));
  ASSERT_EQ(b2, Matrix::packed_diff(b3, b1));
  ASSERT_EQ(b1, Matrix::packed_diff(b3, b2));
  ASSERT_EQ(b0, Matrix::packed_diff(b3, b3));
}

TEST(MatrixTest, PackedMultiplyTest) {
  int8_t a0 = 0x00;
  int8_t a1 = 0x01;
  int8_t a2 = 0x02;
  int8_t a3 = 0x03;
  for (uint8_t i = 0; i < 0xFF; ++i) {
    for (uint8_t j = 0; j < 0xFF; ++j) {
      ASSERT_EQ(Matrix::packed_multiply(i, j), Matrix::packed_multiply(j, i));
    }
  }
  ASSERT_EQ(a0, Matrix::packed_multiply(a0, a0));
  ASSERT_EQ(a0, Matrix::packed_multiply(a0, a1));
  ASSERT_EQ(a0, Matrix::packed_multiply(a0, a2));
  ASSERT_EQ(a0, Matrix::packed_multiply(a0, a3));

  ASSERT_EQ(a0, Matrix::packed_multiply(a1, a0));
  ASSERT_EQ(a1, Matrix::packed_multiply(a1, a1));
  ASSERT_EQ(a2, Matrix::packed_multiply(a1, a2));
  ASSERT_EQ(a3, Matrix::packed_multiply(a1, a3));

  ASSERT_EQ(a0, Matrix::packed_multiply(a2, a0));
  ASSERT_EQ(a2, Matrix::packed_multiply(a2, a1));
  ASSERT_EQ(a0, Matrix::packed_multiply(a2, a2));
  ASSERT_EQ(a2, Matrix::packed_multiply(a2, a3));

  ASSERT_EQ(a0, Matrix::packed_multiply(a3, a0));
  ASSERT_EQ(a3, Matrix::packed_multiply(a3, a1));
  ASSERT_EQ(a2, Matrix::packed_multiply(a3, a2));
  ASSERT_EQ(a1, Matrix::packed_multiply(a3, a3));

  int8_t b0 = a0 | a0 << 2 | a0 << 4 | a0 << 6;
  int8_t b1 = a1 | a1 << 2 | a1 << 4 | a1 << 6;
  int8_t b2 = a2 | a2 << 2 | a2 << 4 | a2 << 6;
  int8_t b3 = a3 | a3 << 2 | a3 << 4 | a3 << 6;

  ASSERT_EQ(b0, Matrix::packed_multiply(b0, b0));
  ASSERT_EQ(b0, Matrix::packed_multiply(b0, b1));
  ASSERT_EQ(b0, Matrix::packed_multiply(b0, b2));
  ASSERT_EQ(b0, Matrix::packed_multiply(b0, b3));

  ASSERT_EQ(b0, Matrix::packed_multiply(b1, b0));
  ASSERT_EQ(b1, Matrix::packed_multiply(b1, b1));
  ASSERT_EQ(b2, Matrix::packed_multiply(b1, b2));
  ASSERT_EQ(b3, Matrix::packed_multiply(b1, b3));

  ASSERT_EQ(b0, Matrix::packed_multiply(b2, b0));
  ASSERT_EQ(b2, Matrix::packed_multiply(b2, b1));
  ASSERT_EQ(b0, Matrix::packed_multiply(b2, b2));
  ASSERT_EQ(b2, Matrix::packed_multiply(b2, b3));

  ASSERT_EQ(b0, Matrix::packed_multiply(b3, b0));
  ASSERT_EQ(b3, Matrix::packed_multiply(b3, b1));
  ASSERT_EQ(b2, Matrix::packed_multiply(b3, b2));
  ASSERT_EQ(b1, Matrix::packed_multiply(b3, b3));
}

TEST(MatrixTest, MatrixAdditionTest) {
  Matrix a({{1, 2, 3}, {4, 5, 6}});
  Matrix b({{3, 2, 1}, {6, 5, 4}});
  Matrix c({{4, 4, 4}, {10, 10, 10}});
  ASSERT_EQ(a + b, c);
}

TEST(MatrixTest, MatrixSubtractionTest) {
  Matrix a({{1, 2, 3}, {4, 5, 6}});
  Matrix b({{3, 2, 1}, {6, 5, 4}});
  Matrix c({{4, 4, 4}, {10, 10, 10}});
  ASSERT_EQ(c - b, a);
}

/**
 * This test compares results of multiplications of 2 random matrices of given size
 * First multiplication is made with trivial algorithm
 * Second multiplication is made with strassen's algorithm
 * Both of results should be equal
 */
TEST(MatrixTest, MatrixMultiplicationTest) {
  size_t sz = 512;

  Matrix f1(sz, sz);
  Matrix f2(sz, sz);
  size_t count = 100;
  for (size_t attempt = 0; attempt < count; ++attempt) {
    std::cout << "\rMultiplication compare " << attempt << "/" << count << std::flush;
    for (size_t i = 0; i < sz; ++i) {
      for (size_t j = 0; j < sz; ++j) {
        f1.set(i, j, rand());
        f2.set(i, j, rand());
      }
    }
    Matrix ff_trivial(Matrix::multiply_trivial(f1, f2));
    Matrix ff_strassen(Matrix::multiply_strassen(f1, f2));
    ASSERT_EQ(ff_trivial, ff_strassen);
  }
  std::cout << "\r";
}
