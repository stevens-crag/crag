/*
 * AAGCryptoTest.cpp
 *
 *  Created on: July 11, 2013
 *      Author: pmorar
 */

#include "gtest/gtest.h"
#include "../include/AAGCrypto.h"

namespace crag {
namespace aag_crypto {


  class Generator {
    public:
      template<typename Rand>
      Generator(const SchemeParameters& params, Rand& rand)
        : k_gen_(params, &rand),
          alice_public_key_(k_gen_.generate_public_key(Participant::Alice)),
          bob_public_key_(k_gen_.generate_public_key(Participant::Bob)),
          alice_private_key_(k_gen_.generate_private_key(alice_public_key_)),
          bob_private_key_(k_gen_.generate_private_key(bob_public_key_)),
          alice_ti_(bob_public_key_, alice_private_key_),
          bob_ti_(alice_public_key_, bob_private_key_) {
//        k_gen_ = make_keys_generator(params, &rand);

//        alice_public_key_ = k_gen_.generate_public_key(Participant::Alice);
//        bob_public_key_ = k_gen_.generate_public_key(Participant::Bob);

//        alice_private_key_ = k_gen_.generate_private_key(alice_public_key_);
//        bob_private_key_ = k_gen_.generate_private_key(bob_public_key_);

//        alice_ti_ = TransmittedInfo(bob_public_key_, alice_private_key_);
//        bob_ti_ = TransmittedInfo(alice_public_key_, bob_private_key_);
      }

      Aut alice_key(CalculationType calc_type) const {
        return k_gen_.make_shared_key(alice_private_key_, bob_ti_, Participant::Alice, calc_type);
      }

      Aut bob_key(CalculationType calc_type) const {
        return k_gen_.make_shared_key(bob_private_key_, alice_ti_, Participant::Bob, calc_type);
      }

    private:
      KeysGenerator k_gen_;
      PublicKey alice_public_key_;
      PublicKey bob_public_key_;

      PrivateKey alice_private_key_;
      PrivateKey bob_private_key_;

      TransmittedInfo alice_ti_;
      TransmittedInfo bob_ti_;
  };

  const int sizes[] = {5};
  const int ITERATIONS_NUM = 5;

  void test_calc_type(CalculationType calc_type) {
    std::default_random_engine rand;
    for (int size: sizes) {
      SchemeParameters params(3, 20, 20, 4, 5, size);

      for (int i = 0; i < ITERATIONS_NUM; ++i) {
        Generator g(params, rand);
        EXPECT_EQ(g.alice_key(calc_type), g.bob_key(calc_type));
      }
    }
  }

  void compare_calc_types(CalculationType calc_type1, CalculationType calc_type2) {
    std::default_random_engine rand;
    for (int size: sizes) {
      SchemeParameters params(3, 20, 20, 4, 5, size);

      for (int i = 0; i < ITERATIONS_NUM; ++i) {
        Generator g(params, rand);
        EXPECT_EQ(g.alice_key(calc_type1), g.alice_key(calc_type2));
      }
    }
  }

  TEST(AAGCryptoTest, SinleReductionCorrectnessTest) {
    test_calc_type(CalculationType::SingleReduction);
  }

  TEST(AAGCryptoTest, IterativeReductionCorrectnessTest) {
    test_calc_type(CalculationType::IterativeReduction);
  }

  TEST(AAGCryptoTest, BlockReductionCorrectnessTest) {
    test_calc_type(CalculationType::BlockReduction);
  }

  TEST(AAGCryptoTest, SinleToIterativeCorrectnessTest) {
    compare_calc_types(CalculationType::SingleReduction, CalculationType::IterativeReduction);
  }

  TEST(AAGCryptoTest, SinleToBlockCorrectnessTest) {
    compare_calc_types(CalculationType::SingleReduction, CalculationType::BlockReduction);
  }

} // namespace fga_crypto
} // namespace crag
