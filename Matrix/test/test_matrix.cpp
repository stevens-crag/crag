#include <gtest/gtest.h>

#include "matrix.h"

namespace crag {
namespace matrix {
namespace {

TEST(Matrix, Ex_01) {
  auto m = Matrix<int>(3);

  EXPECT_EQ(3, m.size1());
  EXPECT_EQ(3, m.size2());
  EXPECT_TRUE(isSquare(m));
  EXPECT_FALSE(isSquare(Matrix<int>(std::make_pair(3, 7))));

  EXPECT_EQ(zero<int>(3), m);

  EXPECT_EQ("[3,3]((0,0,0),(0,0,0),(0,0,0))", m.toString());

  m(0, 0) = 1;

  EXPECT_EQ(1, m(0, 0));
  EXPECT_NE(zero<int>(3), m);
  EXPECT_NE(1, m);

  EXPECT_EQ(1, unit<int>(1));
}

TEST(Matrix, Ex_02) {
  auto m = unit<int>(2);
  EXPECT_EQ(1, m(0, 0));
  EXPECT_EQ(1, m(1, 1));

  m *= 2;
  EXPECT_EQ(2, m(0, 0));
  EXPECT_EQ(2, m(1, 1));

  m += unit<int>(2);
  EXPECT_EQ(3, m(0, 0));
  EXPECT_EQ(3, m(1, 1));

  m -= unit<int>(2);
  EXPECT_EQ(2, m(0, 0));
  EXPECT_EQ(2, m(1, 1));
}

TEST(Matrix, Ex_03) {
  EXPECT_EQ(Matrix<int>(2, {3, 2, 3, 3}), Matrix<int>(2, {1, 2, 3, 4}) + eye<int>({2, -1}));

  EXPECT_THROW({ eye<int>({1}) + eye<int>({1, 2}); }, std::invalid_argument);

  EXPECT_EQ(Matrix<int>(2, {2, 7, 6, 15}), Matrix<int>(2, {1, 2, 3, 4}) * Matrix<int>(2, {2, 1, 0, 3}));

  EXPECT_EQ(Matrix<int>(2, {-1, 1, 3, 1}), Matrix<int>(2, {1, 2, 3, 4}) - Matrix<int>(2, {2, 1, 0, 3}));

  EXPECT_EQ(Matrix<int>(2, {3, 6, 9, 12}), Matrix<int>(2, {1, 2, 3, 4}) * 3);
}


TEST(Matrix, Test_Tr_1) {
  EXPECT_EQ(17, tr(Matrix<int>(2, {2, 7, 6, 15})));
}


TEST(Matrix, Test_CharPoly_1) {
  EXPECT_EQ(std::vector<double>({-40, 4, -10, 1}),  charPoly(Matrix<double>(3, {3, 1, 5, 3, 3, 1, 4, 6, 4})));
  EXPECT_EQ(std::vector<double>({1, -2, 1}),  charPoly(Matrix<double>(2, {2, 1, -1, 0})));
  EXPECT_EQ(std::vector<double>({-2, 5, -4, 1}),  charPoly(Matrix<double>(3, {2, -1, 1, 0, 1, 1, -1, 1, 1})));
  EXPECT_EQ(std::vector<double>({16, -32, 24, -8, 1}),  charPoly(Matrix<double>(4, {1, 1, 0, 0, -1, 3, 0, 0, -6, 8, -1, 1, -16, 22, -9, 5})));
}

} // namespace
} // namespace matrix
} // namespace crag
