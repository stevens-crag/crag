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
    SchemeParameters params(3, 2, 2, 2, 2, 2);
    std::default_random_engine rnd;
    KeysGenerator alice(SchemeParameters::canonical(), rnd);
  }

} // namespace fga_crypto
} // namespace crag
