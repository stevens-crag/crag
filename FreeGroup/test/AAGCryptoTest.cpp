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

  template<typename Rand>
  std::pair<Aut, Aut> simulate_key_exchange(const SchemeParameters& params, Rand& rand) {
    auto k_gen = make_keys_generator(params, &rand);

    auto a_pub = k_gen.generate_public_key(Alice);
    auto b_pub = k_gen.generate_public_key(Bob);

    auto a_priv = k_gen.generate_private_key(a_pub);
    auto b_priv = k_gen.generate_private_key(b_pub);

    TransmittedInfo a_ti(b_pub, a_priv);
    TransmittedInfo b_ti(a_pub, b_priv);

    auto a_shared = k_gen.make_shared_key(a_priv, b_ti, Alice);
    auto b_shared = k_gen.make_shared_key(b_priv, a_ti, Bob);

    std::cout << "as: vn=" << slp_vertices_num(a_shared) <<
                 ", bs: vn=" << slp_vertices_num(b_shared) << std::endl;

    auto a_fr = a_shared.free_reduction();
    auto b_fr = b_shared.free_reduction();
    return std::make_pair(a_fr, b_fr);
  }


  TEST(AAGCryptoTest, CorrectnessTest) {
    std::default_random_engine rand;
    for (int size = 0; size < 5; ++size) {
      SchemeParameters params(3, size, size, size, size, size);

      for (int i = 0; i < 5; ++i) {
        auto keys = simulate_key_exchange(params, rand);

        EXPECT_EQ(keys.first, keys.second);
      }
    }
  }

  TEST(AAGCryptoTest, SpeedTest) {
    std::default_random_engine rand;
    SchemeParameters params(3, 20, 20, 5, 8, 4);

    auto keys = simulate_key_exchange(params, rand);

    EXPECT_EQ(keys.first, keys.second);
  }

} // namespace fga_crypto
} // namespace crag
