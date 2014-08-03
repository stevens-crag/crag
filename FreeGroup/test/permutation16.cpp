/*
 * permutation16.cpp
 *
 */

#include <vector>
#include <algorithm>
#include <random>
#include <unordered_map>

#include "gtest/gtest.h"
#include "permutation16.h"
#include "Permutation.h"

namespace crag {
namespace {

TEST(Permutation16, TrivialPermuation) {
  Permutation16 permutation;
  for (size_t i = 0; i < 16; ++i) {
    EXPECT_EQ(i, permutation[i]);
  }
}

TEST(Permutation16, SetElement) {
  Permutation16 permutation;
  permutation.set_image(7, 9);
  for (size_t i = 0; i < 16; ++i) {
    if (i != 7) {
      EXPECT_EQ(i, permutation[i]);
    } else {
      EXPECT_EQ(9, permutation[i]);
    }
  }

  Permutation16 pre_defined = {11, 4, 5, 13, 15, 0, 12, 8, 3, 1, 6, 14, 9, 7, 2, 10};
  size_t element = 0;
  for (auto expected : {11, 4, 5, 13, 15, 0, 12, 8, 3, 1, 6, 14, 9, 7, 2, 10}) {
    EXPECT_EQ(expected, pre_defined[element]) << element;
    ++element;
  }

  Permutation16 pre_defined_short = {11, 4, 5, 13, 0, 12, 8, 3, 1, 6, 14, 9, 7, 2, 10};
  element = 0;
  for (auto expected : {11, 4, 5, 13, 0, 12, 8, 3, 1, 6, 14, 9, 7, 2, 10, 15}) {
    EXPECT_EQ(expected, pre_defined_short[element]) << element;
    ++element;
  }

}

TEST(Permutation16, Swap) {
  for (size_t i = 0; i < 16; ++i) {
    for (size_t j = i; j < 16; ++j) {
      Permutation16 permutation;
      permutation.swap_images(i, j);
      for (size_t k = 0; k < 16; ++k) {
        if (k != i && k != j) {
          EXPECT_EQ(k, permutation[k]);
        } else  if (k == i) {
          EXPECT_EQ(j, permutation[k]);
        } else  if (k == j) {
          EXPECT_EQ(i, permutation[k]);
        }
      }
    }
  }
}

TEST(Permutation16, ComposeAndInverseCheck) {
  for (size_t first_size = 2; first_size <= 6; ++first_size) {
    std::vector<uint64_t> first_permutation(first_size);
    for (size_t i = 0; i < first_size; ++i) {
      first_permutation[i] = i;
    }
    for (size_t second_size = 2; second_size <= 6; ++second_size) {
      std::vector<uint64_t> second_permutation(second_size);
      for (size_t i = 0; i < second_size; ++i) {
        second_permutation[i] = i;
      }
      do {
        Permutation first(first_permutation.begin(), first_permutation.end());
        Permutation16 first16(first_permutation.begin(), first_permutation.end());

//        std::stringstream first_string;
//        for (size_t i = 0; i < first_permutation.size(); ++i) {
//          first_string << first_permutation[i] << ',';
//        }

        for (size_t i = 0; i < first_size; ++i) {
          ASSERT_EQ(first[i], first16[i]);
        }

        Permutation first_inverse = first.inverse();
        Permutation16 first16_inverse = first16.inverse();
        for (size_t i = 0; i < first_size; ++i) {
          ASSERT_EQ(first_inverse[i], first16_inverse[i]);
        }

        do {
          Permutation second(second_permutation.begin(), second_permutation.end());
          Permutation16 second16(second_permutation.begin(), second_permutation.end());
          std::stringstream second_string;
//          for (size_t i = 0; i < second_permutation.size(); ++i) {
//            second_string << second_permutation[i] << ',';
//          }

          auto composition = first * second;
          auto composition16 = first16 * second16;

//          std::stringstream composition_string;
//          for (size_t i = 0; i < composition.size(); ++i) {
//            composition_string << composition[i] << ',';
//          }
//
//          std::stringstream composition16_string;
//          for (size_t i = 0; i < composition.size(); ++i) {
//            composition16_string << composition16[i] << ',';
//          }

          ASSERT_EQ(static_cast<size_t>(composition.size()), composition16.size());

          for (size_t i = 0; i < composition16.size(); ++i) {
            ASSERT_EQ(composition[i], composition16[i]);// << first_string.str() << " * " << second_string.str() << std::endl
                //<< composition_string.str() << " or " << composition16_string.str();
          }

          for (size_t i = composition16.size(); i < 16; ++i) {
            ASSERT_EQ(i, composition16[i]);// << composition16.permutation_;
          }
        } while (std::next_permutation(second_permutation.begin(), second_permutation.end()));
      } while (std::next_permutation(first_permutation.begin(), first_permutation.end()));
    }
  }
}

constexpr uint64_t factorial(uint64_t n)
{
    return n > 0? factorial(n-1) * n : 1;
}

TEST(Permutation16, SmallRandom) {
  const size_t REPEAT = 10000000ul;
  const size_t PERMUTATION_SIZE = 6;
  mt19937 generator(16);

  std::unordered_map<Permutation16, size_t> hash_count;
  hash_count.reserve(factorial(PERMUTATION_SIZE));

  for (size_t count = 0; count < REPEAT; ++count) {
    hash_count[Permutation16::random(PERMUTATION_SIZE - 1, generator)]++;
  }

  ASSERT_EQ(factorial(PERMUTATION_SIZE), hash_count.size());

  double ideal_count = static_cast<double>(REPEAT) / factorial(PERMUTATION_SIZE);
  for (auto count : hash_count) {
    ASSERT_GE(count.second, ideal_count * 0.95);
    ASSERT_LE(count.second, ideal_count * 1.05);
  }
}

TEST(Permutation16, LargeRandom) {
  const size_t REPEAT = 10000000ul;
  const size_t PERMUTATION_SIZE = 16;
  mt19937 generator(16);

  size_t image_count_0[PERMUTATION_SIZE] = {0};
  size_t image_count_15[PERMUTATION_SIZE] = {0};

  for (size_t count = 0; count < REPEAT; ++count) {
    auto permutation = Permutation16::random(PERMUTATION_SIZE - 1, generator);
    ++image_count_0[permutation[0]];
    ++image_count_15[permutation[PERMUTATION_SIZE - 1]];
  }

  double ideal_count = static_cast<double>(REPEAT) / PERMUTATION_SIZE;
  for (auto count : image_count_0) {
    ASSERT_GE(count, ideal_count * 0.95);
    ASSERT_LE(count, ideal_count * 1.05);
  }
  for (auto count : image_count_15) {
    ASSERT_GE(count, ideal_count * 0.95);
    ASSERT_LE(count, ideal_count * 1.05);
  }
}


}
}

