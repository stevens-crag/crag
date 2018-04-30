#include "gtest/gtest.h"

#include "FiniteField.h"

namespace crag {
namespace finitefield {
namespace {

TEST(FiniteField, PrimeFieldTest1) {
  using ZZ5 = ZZ<5>;

  const ZZ5 zero(0);
  const ZZ5 one(1);
  const ZZ5 two(2);
  const ZZ5 three(3);
  const ZZ5 four(4);

  EXPECT_NE(zero, one);
  EXPECT_EQ(one, one);
  EXPECT_EQ(one, two * three);
  EXPECT_EQ(two, three + four);
  EXPECT_EQ(four, three - four);
  EXPECT_EQ(two, one / three);
  EXPECT_EQ(two, pwr(three, -1));
  EXPECT_EQ(four, pwr(three, 2));
  EXPECT_EQ(four, pwr(three, 10));
  EXPECT_EQ(two, pwr(three, 11));
}

TEST(FiniteField, FiniteFieldTest1) {
  using GF256 = FieldElement<IdealGeneratedByPolynomial<ZZ<2>, 1, 1, 0, 1, 1, 0, 0, 0, 1>>;

  const GF256 zero(0);
  const GF256 one(1);
  const GF256 a({1, 1, 0, 0, 1, 0, 1});
  const GF256 b({0, 1, 0, 1, 0, 0, 1, 1});
  const GF256 c({1, 1});
  const GF256 d({1, 0, 1});

  EXPECT_EQ(zero, one + one);
  EXPECT_EQ(one, a * b);
  EXPECT_EQ(b, pwr(a, -1));
  EXPECT_EQ(d, pwr(c, 2));
  EXPECT_EQ(b, one / a);
  EXPECT_THROW(a / zero, std::logic_error);
}

TEST(FiniteField, RandomTest1) {
  using Z5 = FieldElement<IdealGeneratedByInteger<5>>;

  std::mt19937 g(1233);

  EXPECT_EQ(3, Z5::random(g));
}

TEST(FiniteField, RandomTest2) {
  using Z2 = FieldElement<IdealGeneratedByInteger<2>>;
  using GF256 = FieldElement<IdealGeneratedByPolynomial<Z2, 1, 1, 0, 1, 1, 0, 0, 0, 1>>;

  std::mt19937 g(1233);

  EXPECT_EQ(GF256({1, 1, 1, 0, 1, 0, 1, 1, 1}), GF256::random(g));
}
} // namespace
} // namespace finitefield
} // namespace crag
