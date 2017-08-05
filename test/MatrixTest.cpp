#include <stdexcept>

#include <gtest/gtest.h>

#include "matrix_strassen.h"

class MatrixTest : public ::testing::Test
{
protected:
  void SetUp()
  {
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

TEST(MatrixTest, MatrixMultiplicationTest) {
  Matrix a({{1, 2, 3}, {4, 5, 6}});
  Matrix b({{1, 4}, {2, 5}, {3, 6}});
  Matrix c({{14, 32}, {32, 77}});
  Matrix d(100, 100);
  Matrix e(100, 100);
  for (size_t i = 0; i < 100; ++i) {
    for (size_t j = 0; j < 100; ++j) {
      d.set(i, j, 1);
      e.set(i, j, 100);
    }
  }
  Matrix f(1024, 1024);
  for (size_t i = 0; i < 1024; ++i) {
    f.set(i, i, 1);
  }
  ASSERT_EQ(a*b, c);
  ASSERT_EQ(d*d, e);
  ASSERT_EQ(f*f, f);
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
