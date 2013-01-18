#include "gtest/gtest.h"
#include "SLPSet.h"

//! The fixture for testing class SLPSet.
class SLPSetTest : public ::testing::Test {
 protected:
//TODO fill methods
  FooTest() {
    // You can do set-up work for each test here.
  }

  virtual ~FooTest() {
    // You can do clean-up work that doesn't throw exceptions here.
  }

  // If the constructor and destructor are not enough for setting up
  // and cleaning up each test, you can define the following methods:

  virtual void SetUp() {
    // Code here will be called immediately after the constructor (right
    // before each test).
  }

  virtual void TearDown() {
    // Code here will be called immediately after each test (right
    // before the destructor).
  }

  // Objects declared here can be used by all tests in the test case for Foo.
};

TEST(SLPSetTest, DefaultSet) {
  SLPSet s1(5);
  SLPSet s2(5);
  SLPSet s3(6);
  EXPECT_EQ(s1, s2);
  EXPECT_NE(s1, s3);
  SLPSet&& freeReductions = s2.free_reduction();

}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

