/*
 * SLPSet_main.cpp
 *
 *  Created on: Jan 21, 2013
 *      Author: pmorar
 */

#include "gtest/gtest.h"
#include "SLPSet.h"

namespace crag {

/**
 * Wrap for SLPProducedWordIterator, which returns printable string for the vertex.
 */
class Iterator: public std::iterator<std::forward_iterator_tag, //iterator_category
const std::string                 //value_type
>  {
	SLPProducedWord::const_iterator it_;
public:
	Iterator(SLPProducedWord::const_iterator it) {
		it_ = it;
	}

	Iterator& operator++() {
		it_++;
		return *this;
	}

	Iterator operator++(int) { //!< Postincrement
		Iterator copy(*this);
		++(*this);
		return copy;
	}

	std::string operator*() const { //!< "Dereference" current symbol
		auto v = *it_;
		std::stringstream s(v.is_negative() ? "-" : "+");
		s << v.terminal_symbol();
		return s.str();
	}

	bool operator==(const Iterator& other) const {
		return it_ == other.it_;
	}

	bool operator!=(const Iterator& other) const {
		return !(*this == other);
	}
};

}

namespace std {

::std::ostream& operator<<(::std::ostream& os, const ::crag::SLPVertex& vertex)
{
    os << "vertex " << "(" << vertex.length().get_str() << ", " << vertex.height()
        << ", " << vertex.terminal_symbol() << ") - " << ::testing::PrintToString(vertex) << ::std::endl;
    os << "w[";
    crag::SLPProducedWord w(vertex);
	copy(crag::Iterator(w.begin()), crag::Iterator(w.end()), ostream_iterator<string>(os));
	os << "]" << endl;
	return os;
}
}

namespace crag {

class SLPSetBasicTest : public ::testing::Test {
    protected:
	SLPSetBasicTest()
        : v(SLPVertex::concatenate(
        		SLPVertex::terminal_vertex(1),
        		SLPVertex::concatenate(
        				SLPVertex::terminal_vertex(2),
        				SLPVertex::terminal_vertex(1)))),
       u(SLPVertex::terminal_vertex(5)),
       slp(v),
       slp2({u, v})
      {
      }

    SLPVertex v, u;
  	SLPSet slp, slp2;


  };

TEST_F(SLPSetBasicTest, SimpleConstructor) {
for (int i = 0, n = 0; i < 10; ++i, n += i) {
  SLPSet slp(n);
  ASSERT_EQ(slp.roots_num(), n);
  for (int j = 0; j < n; ++j) {
	SLPProducedWord word = slp.produced_word(j);
	EXPECT_EQ(word.size(), 1);
	const SLPVertex& v = word[0];
	EXPECT_TRUE(v.is_terminal());
	EXPECT_EQ(v.terminal_symbol(), j + 1);
  }
}
}

TEST_F(SLPSetBasicTest, SingleRootConstructor) {
	EXPECT_EQ(slp.roots_num(), 1);
	EXPECT_EQ(slp.root(0), v);
}

TEST_F(SLPSetBasicTest, RangeConstructor) {
	std::vector<SLPVertex> range(10, v);
	SLPSet slp1(range.begin(), range.end());
	EXPECT_EQ(slp1.roots_num(), 10);
	for (int i = 0; i < 10; ++i)
		EXPECT_EQ(slp1.root(i), v) << "i = " << i;
}

TEST_F(SLPSetBasicTest, InitializerListConstructor) {
	EXPECT_EQ(slp2.roots_num(), 2);
	SLPSet slp1 = {v, v, v, v, v};
	EXPECT_EQ(slp1.roots_num(), 5);
	for (int i = 0; i < 5; ++i)
		EXPECT_EQ(slp1.root(i), v);
}

TEST_F(SLPSetBasicTest, SimpleCopyConstructor) {
	SLPSet slp1(slp);
	EXPECT_EQ(slp1.roots_num(), 1);
	EXPECT_EQ(slp.root(0), v);

	SLPSet slp21(slp2);
	EXPECT_EQ(slp21.roots_num(), 2);
	EXPECT_EQ(slp21.root(0), u);
	EXPECT_EQ(slp21.root(1), v) << "slp21.root(0) = " << slp21.root(0)
					<< "v = " << v;
}

TEST_F(SLPSetBasicTest, SimpleMoveConstructor) {
	SLPSet slp1(std::move(slp));
	EXPECT_EQ(slp1.roots_num(), 1);
	EXPECT_EQ(slp.roots_num(), 0);
	EXPECT_EQ(slp1.root(0), v);

	SLPSet slp21(std::move(slp2));
	EXPECT_EQ(slp21.roots_num(), 2);
	EXPECT_EQ(slp2.roots_num(), 0);
	EXPECT_EQ(slp21.root(0), u);
	EXPECT_EQ(slp21.root(1), v);
}



TEST_F(SLPSetBasicTest, EqualityOperator) {
	EXPECT_EQ(slp, slp) << "SLPSet does not equal itself";
	EXPECT_EQ(SLPSet(5), SLPSet(5)) << "SLPSets consisting of the same terminals must be equal";
	EXPECT_NE(SLPSet(5), SLPSet(4)) << "SLPSets consisting of different number of terminals must be non-euqal";
	EXPECT_NE(SLPSet(1), slp) << "One root terminal equals one root non trivial";


	EXPECT_NE(SLPSet({v, u}),
			SLPSet({u, v})) << "Order of vertices should matter";

	EXPECT_EQ(SLPSet({v, v, u, v, u}), SLPSet({v, v, u, v, u}))
					<< "SLPSets with the same content must be equal";
}

TEST_F(SLPSetBasicTest, RootInvertionTest) {
	auto new_v = slp.invert_root(0).root(0);
	EXPECT_EQ(new_v, v.negate());
}

TEST_F(SLPSetBasicTest, RootReplacementTest) {
	auto new_v = slp.replace_root(0, u).root(0);
	EXPECT_EQ(new_v, u);
}

 
}
