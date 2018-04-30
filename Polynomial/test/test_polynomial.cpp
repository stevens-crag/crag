#include <gtest/gtest.h>

#include "polynomial.h"

namespace crag {
namespace polynomials {
namespace {

template <int p>
using GF = finitefield::FieldElement<finitefield::IdealGeneratedByInteger<p>>;
using GF5 = GF<5>;

TEST(Monomial, BasicChecks) {
  EXPECT_THROW({ IntStandardMonomial(IntStandardMonomial::term_t({})); }, std::invalid_argument);
  EXPECT_THROW({ IntStandardMonomial({0}); }, std::invalid_argument);

  EXPECT_TRUE(IntStandardMonomial({0, 0}).isUnit());
  EXPECT_TRUE(IntStandardMonomial({0, 0}).isNumber());
  EXPECT_TRUE(1 == IntStandardMonomial({0, 0}));

  EXPECT_FALSE(IntStandardMonomial({0, 0}).isZero());
  EXPECT_TRUE(IntStandardMonomial({0, 0}) != 0);

  EXPECT_TRUE(IntStandardMonomial(0, {1, 1}).isZero());
  EXPECT_TRUE(IntStandardMonomial(0, {1, 1}).isNumber());
  EXPECT_TRUE(0 == IntStandardMonomial(0, {1, 1}));
  EXPECT_FALSE(0 != IntStandardMonomial(0, {1, 1}));

  const auto m = IntStandardMonomial({1, 2, 3});

  EXPECT_EQ(3, m.dimension());
  EXPECT_EQ(IntStandardMonomial::term_t({1, 2, 3}), m);
  EXPECT_EQ(IntStandardMonomial::term_t({1, 2, 3}), m.term());
  EXPECT_EQ(1, m.coef());

  EXPECT_FALSE(m.isZero());
  EXPECT_FALSE(m.isUnit());

  EXPECT_EQ(2, IntStandardMonomial(2, {0, 0, 0}));
  EXPECT_FALSE(IntStandardMonomial(2, {0, 3, 0}).isNumber());

  // zeros are equal
  EXPECT_EQ(IntStandardMonomial(0, {1, 2}), IntStandardMonomial(0, {3, 5}));

  EXPECT_THROW({ IntStandardMonomial(0, {1, 2}) == IntStandardMonomial(0, {3, 5, 1}); }, std::invalid_argument);
}

TEST(Monomial, Evaluation_1) {
  const auto m = StandardMonomial<GF5>(GF5(2), {2, 3, 1});
  const auto value = evaluate(m, {GF5(2), GF5(3), GF5(4)});

  EXPECT_EQ(GF5(4), value);

  EXPECT_EQ(GF5(0), evaluate(StandardMonomial<GF5>(GF5(0), {2, 3}), {GF5(2), GF5(3)}));
}

TEST(Monomial, Evaluation_2) {
  const auto m = LaurentMonomial<GF5>(GF5(2), {2, -3, 1});
  EXPECT_THROW({ evaluate(m, {GF5(2), GF5(0), GF5(4)}); }, std::logic_error);
}

TEST(Monomial, Operations) {
  const auto m = IntStandardMonomial({1, 2, 3});

  EXPECT_EQ(IntStandardMonomial(3, {1, 2, 3}), m * 3);
  EXPECT_EQ(IntStandardMonomial(5, {1, 2, 3}), 5 * m);

  EXPECT_EQ(IntStandardMonomial({3, 5, 4}), m * IntStandardMonomial::term_t({2, 3, 1}));
  EXPECT_EQ(IntStandardMonomial(4, {3, 5, 4}), m * IntStandardMonomial(4, {2, 3, 1}));

  EXPECT_EQ(IntStandardMonomial({3, 5, 4}), IntStandardMonomial::term_t({2, 3, 1}) * m);
}

TEST(LaurentMonomial, Operations) {
  auto m = IntLaurentMonomial({1, 2, 3});

  m *= 2;
  EXPECT_EQ(2, m.coef());

  m *= IntLaurentMonomial::term_t({1, 0, 3});
  EXPECT_EQ(IntLaurentMonomial(2, {2, 2, 6}), m);

  m /= IntLaurentMonomial::term_t({3, 5, 3});
  EXPECT_EQ(IntLaurentMonomial(2, {-1, -3, 3}), m);

  EXPECT_EQ(IntLaurentMonomial(2, {-2, -4, 0}), m / IntLaurentMonomial::term_t({1, 1, 3}));

  EXPECT_THROW({ m* IntLaurentMonomial({1, 2}); }, std::invalid_argument);
}

TEST(Polynomial, BasicChecks) {
  const auto p = IntPolynomial(2);

  EXPECT_TRUE(p.isZero());
  EXPECT_FALSE(p.isUnit());

  EXPECT_EQ(0, p.size());
  EXPECT_EQ(2, p.dimension());

  EXPECT_EQ(p.end(), p.begin());

  EXPECT_EQ("0", p.toString());
  EXPECT_EQ("1", (p + 1).toString());
  EXPECT_EQ("1 + x_1 * x_2", (p + IntStandardMonomial({0, 0}) + IntStandardMonomial({1, 1})).toString());

  EXPECT_THROW({ IntPolynomial(1); }, std::invalid_argument);
}

TEST(Polynomial, Sum) {
  auto p = IntPolynomial(2);

  EXPECT_EQ("0", (p + IntStandardMonomial(0, {1, 1})).toString());
  EXPECT_EQ("x_1 * x_2", (p + IntStandardMonomial({1, 1})).toString());
  EXPECT_EQ("x * y", (p + IntStandardMonomial({1, 1})).toString({"x", "y"}));
  EXPECT_EQ("2 * x_1 * x_2^2", (IntStandardMonomial(2, {1, 2}) + p).toString());

  EXPECT_EQ("x_1 * x_2", (p + IntPolynomial::term_t({1, 1})).toString());
  EXPECT_EQ("x_1 * x_2^2", (IntPolynomial::term_t({1, 2}) + p).toString());

  EXPECT_EQ("2", (p + 2).toString());
  EXPECT_EQ("3", (3 + p).toString());
  EXPECT_EQ("0", (0 + p).toString());

  p += IntStandardMonomial({1, 1});
  EXPECT_EQ("0", (p + IntStandardMonomial(-1, {1, 1})).toString());
}

TEST(Polynomial, Difference) {
  auto p = IntPolynomial(2);

  p += IntStandardMonomial(2, {1, 2});

  EXPECT_EQ("2 * x_1 * x_2^2", p.toString());
  EXPECT_EQ("-2 * x_1 * x_2^2", (-p).toString());

  EXPECT_EQ(p.toString(), (p - IntStandardMonomial(0, {1, 1})).toString());
  EXPECT_EQ("-1 * x_1 * x_2 + 2 * x_1 * x_2^2", (p - IntStandardMonomial({1, 1})).toString());
  EXPECT_EQ("0", (IntStandardMonomial(2, {1, 2}) - p).toString());

  EXPECT_EQ("-1 * x_1 * x_2 + 2 * x_1 * x_2^2", (p - IntPolynomial::term_t({1, 1})).toString());
  EXPECT_EQ("-1 * x_1 * x_2^2", (IntPolynomial::term_t({1, 2}) - p).toString());

  EXPECT_EQ("-2 + 2 * x_1 * x_2^2", (p - 2).toString());
  EXPECT_EQ("3 + -2 * x_1 * x_2^2", (3 - p).toString());
  EXPECT_EQ("2 * x_1 * x_2^2", (p - 0).toString());
  EXPECT_EQ("-2 * x_1 * x_2^2", (0 - p).toString());
}

TEST(Polynomial, Product_1) {
  auto p = IntPolynomial(2);

  p += IntStandardMonomial(2, {1, 2});

  EXPECT_EQ("0", (p * IntStandardMonomial(0, {1, 1})).toString());
  EXPECT_EQ("2 * x_1^2 * x_2^3", (p * IntStandardMonomial({1, 1})).toString());
  EXPECT_EQ("4 * x_1^2 * x_2^4", (IntStandardMonomial(2, {1, 2}) * p).toString());

  EXPECT_EQ("2 * x_1^2 * x_2^3", (p * IntPolynomial::term_t({1, 1})).toString());
  EXPECT_EQ("2 * x_1^2 * x_2^4", (IntPolynomial::term_t({1, 2}) * p).toString());

  EXPECT_EQ("4 * x_1 * x_2^2", (p * 2).toString());
  EXPECT_EQ("6 * x_1 * x_2^2", (3 * p).toString());
  EXPECT_EQ("0", (p * 0).toString());
}

TEST(Polynomial, Product_2) {
  auto p = IntPolynomial(2);
  p = p + IntStandardMonomial({1, 0}) + IntStandardMonomial({0, 1});

  EXPECT_EQ("x_2 + x_1", p.toString());

  EXPECT_EQ("x_2^2 + 2 * x_1 * x_2 + x_1^2", (p * p).toString());
}

TEST(LaurentPolynomial, DivisionByTerm) {
  auto p = IntLaurentPolynomial(2);

  p = p + IntLaurentMonomial(2, {1, 2}) + IntLaurentMonomial(3, {2, -2});

  EXPECT_FALSE(p.isZero());
  EXPECT_FALSE(p.isUnit());

  EXPECT_EQ("2 * x_1 * x_2^2 + 3 * x_1^2 * x_2^-2", p.toString());

  EXPECT_EQ("2 * x_2 + 3 * x_1 * x_2^-3", (p / IntLaurentPolynomial::term_t({1, 1})).toString());
}

TEST(Polynomial, Comparison) {
  auto p = IntPolynomial(1, IntPolynomial::term_t({1, 0}));

  EXPECT_FALSE(p.isZero());
  EXPECT_FALSE(p.isUnit());
  EXPECT_FALSE(p.isNumber());

  EXPECT_THROW({ p == IntStandardMonomial({1, 0, 9}); }, std::invalid_argument);

  EXPECT_EQ(p, IntStandardMonomial({1, 0}));
  EXPECT_EQ(IntStandardMonomial({1, 0}), p);
  EXPECT_NE(p, IntStandardMonomial({1, 1}));
  EXPECT_NE(IntStandardMonomial({1, 1}), p);

  EXPECT_EQ(p, IntPolynomial::term_t({1, 0}));
  EXPECT_EQ(IntPolynomial::term_t({1, 0}), p);
  EXPECT_NE(p, IntPolynomial::term_t({1, 1}));
  EXPECT_NE(IntPolynomial::term_t({1, 1}), p);

  EXPECT_NE(p, 2);
  EXPECT_NE(2, p);

  // compare zeros
  EXPECT_EQ(IntPolynomial(2), IntStandardMonomial(0, {1, 1}));

  EXPECT_THROW({ IntPolynomial(2) == IntPolynomial(3); }, std::invalid_argument);

  EXPECT_TRUE(IntPolynomial(2).isZero());
  EXPECT_TRUE(IntPolynomial(2).isNumber());
}

TEST(Polynomial, Evaluation_1) {
  const auto p = Polynomial<GF5>(2);
  EXPECT_EQ(GF5(0), evaluate(p, {GF5(2), GF5(3)}));

  const auto q = p + StandardMonomial<GF5>({2, 3}) + StandardMonomial<GF5>(3, {1, 2});
  EXPECT_EQ(GF5(3), evaluate(q, {GF5(3), GF5(2)}));
}

TEST(Polynomial, Evaluation_2) {
  const auto p = LaurentPolynomial<GF5>(2) + LaurentMonomial<GF5>({-2, 3});
  EXPECT_THROW({ evaluate(p, {GF5(0), GF5(1)}); }, std::logic_error);
}
} // namespace
} // namespace polynomials
} // namespace crag