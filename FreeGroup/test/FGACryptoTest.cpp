/*
 * FGACrypto.h
 *
 *  Created on: May 2, 2013
 *      Author: pmorar
 */

#include "gtest/gtest.h"
#include "../include/FGACrypto.h"

namespace crag {
namespace fga_crypto {
//  TEST(FGACryptoTest, Basic) {
////    gmp_pool_setup();
//    SchemeParameters params(3, 5, 5, 3, 3, 1);
//    std::default_random_engine rnd(2235);
//    KeysGenerator alice(params, rnd);

//    auto trivial_shared_key = alice.make_shared_key(alice.public_keys());

//    EXPECT_EQ(EndomorphismSLP<int>::identity(), trivial_shared_key());
//  }

  class MockParticipant {
    public:
      MockParticipant(const AutDescription& private_key, const AutDescription& other_private_key)
        : priv_key_(private_key),
          other_priv_key_(other_private_key) {}

      PublicKeys process_incoming_public_keys(const PublicKeys& incoming_public_keys) {
        PublicKeys result;
        result.s = conjugate_all(incoming_public_keys.s, priv_key_);
        result.r = conjugate_all(incoming_public_keys.r, priv_key_);
        return result;
      }

      AutDescription make_shared_key(bool isAlice = true) {
        if (isAlice) {//we are Alice
          return  priv_key_ * other_priv_key_ * priv_key_.inverse_description() * other_priv_key_.inverse_description();//a*b*a^{-1}*b^{-1}
        } else {//we are Bob
          return other_priv_key_ * priv_key_ * other_priv_key_.inverse_description() * priv_key_.inverse_description();
        }
      }

    private:
      AutDescription priv_key_;
      AutDescription other_priv_key_;
  };

  TEST(FGACryptoTest, CheckMock) {
//    gmp_pool_setup();
    std::default_random_engine rnd(2235);
    UniformAutomorphismSLPGenerator<int> aut_random(3, &rnd);
    AutDescription a_private_key(3, aut_random);
    AutDescription b_private_key(3, aut_random);

    MockParticipant alice(a_private_key, b_private_key);
    MockParticipant bob(b_private_key, a_private_key);

    auto a_shared = alice.make_shared_key();
    auto b_shared = bob.make_shared_key(false);

    EXPECT_EQ(a_shared(), b_shared());
  }

  TEST(FGACryptoTest, MockShared) {
//    gmp_pool_setup();
    SchemeParameters params(3, 5, 5, 3, 3, 2);
    std::default_random_engine rnd(2235);
    KeysGenerator real(params, rnd);
    for (int i = 0; i < 10; ++i) {
      UniformAutomorphismSLPGenerator<int> aut_random(3, &rnd);
      AutDescription mock_private_key(3, aut_random);

      MockParticipant mock(mock_private_key, real.private_key());

      auto processed_pub_key = mock.process_incoming_public_keys(real.public_keys());

      auto alice_shared_key = real.make_shared_key(processed_pub_key);
      auto bob_shared_key = mock.make_shared_key(false);

      EXPECT_EQ(alice_shared_key().free_reduction(), bob_shared_key().free_reduction());

      alice_shared_key = mock.make_shared_key(true);
      bob_shared_key = real.make_shared_key(processed_pub_key, false);

      EXPECT_EQ(alice_shared_key().free_reduction(), bob_shared_key().free_reduction());
    }
  }

  TEST(FGACryptoTest, SharedKeys) {
//    gmp_pool_setup();
    SchemeParameters params(3, 2, 2, 2, 2, 1);
    std::default_random_engine rnd(22345);
    for (int i = 0; i < 5; ++i) {

//      std::cout << "start" << i << std::endl;
      KeysGenerator alice(params, rnd);
      KeysGenerator bob(params, rnd);

      auto a_pk = alice.public_keys();
      auto b_pk = bob.public_keys();

      auto a_processed = bob.process_incoming_public_keys(a_pk);
      auto b_processed = alice.process_incoming_public_keys(b_pk);

      auto a_shared_key = alice.make_shared_key(a_processed, true);
      auto b_shared_key = bob.make_shared_key(b_processed, false);

//      std::cout << "comparing keys" << std::endl;
      EXPECT_EQ(a_shared_key().free_reduction(), b_shared_key().free_reduction());
//      std::cout << "finish" << i << std::endl;

      //basic speed for 3, 5, 4, 2, 2, 1 is 6.7s
      //better 2.9s
    }
  }

} // namespace fga_crypto
} // namespace crag
