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
  TEST(FGACryptoTest, Basic) {
//    gmp_pool_setup();
    SchemeParameters params(3, 5, 5, 3, 3, 1);
    std::default_random_engine rnd(2235);
    KeysGenerator alice(params, rnd);

    auto trivial_shared_key = alice.make_shared_key(alice.public_keys());

    EXPECT_EQ(EndomorphismSLP<int>::identity(), trivial_shared_key());
  }

//  TEST(FGACryptoTest, SharedKeys) {
////    gmp_pool_setup();
//    SchemeParameters params(3, 1, 1, 1, 1, 0);
//    std::default_random_engine rnd(22345);
//    KeysGenerator alice(params, rnd);
//    KeysGenerator bob(params, rnd);

//    auto a_pk = alice.public_keys();
//    auto b_pk = bob.public_keys();

//    auto a_processed = bob.process_incoming_public_keys(a_pk);
//    auto b_processed = alice.process_incoming_public_keys(b_pk);

//    auto a_shared_key = alice.make_shared_key(a_processed);
//    auto b_shared_key = alice.make_shared_key(b_processed, false);

//    EXPECT_EQ(a_shared_key().free_reduction(), b_shared_key().free_reduction());
//  }

} // namespace fga_crypto
} // namespace crag
