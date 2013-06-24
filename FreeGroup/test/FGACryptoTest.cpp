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
      MockParticipant(const AutomorphismDescription& private_key, const AutomorphismDescription& other_private_key)
        : priv_key_(private_key),
          other_priv_key_(other_private_key) {}

      PublicKeys process_incoming_public_keys(const PublicKeys& incoming_public_keys) {
        PublicKeys result;
        result.s = conjugate_all(incoming_public_keys.s, priv_key_);
        result.r = conjugate_all(incoming_public_keys.r, priv_key_);
        return result;
      }

      AutomorphismDescription make_shared_key(bool isAlice = true) {
        if (isAlice) {//we are Alice
          return  priv_key_ * other_priv_key_ * priv_key_.inverse_description() * other_priv_key_.inverse_description();//a*b*a^{-1}*b^{-1}
        } else {//we are Bob
          return other_priv_key_ * priv_key_ * other_priv_key_.inverse_description() * priv_key_.inverse_description();
        }
      }

    private:
      AutomorphismDescription priv_key_;
      AutomorphismDescription other_priv_key_;
  };

//  TEST(FGACryptoTest, CheckMock) {
////    gmp_pool_setup();
//    std::default_random_engine rnd(2235);
//    UniformAutomorphismSLPGenerator<int> aut_random(3, &rnd);
//    AutomorphismDescription a_private_key(3, aut_random);
//    AutomorphismDescription b_private_key(3, aut_random);

//    MockParticipant alice(a_private_key, b_private_key);
//    MockParticipant bob(b_private_key, a_private_key);

//    auto a_shared = alice.make_shared_key();
//    auto b_shared = bob.make_shared_key(false);

//    EXPECT_EQ(a_shared(), b_shared());
//  }

//  TEST(FGACryptoTest, MockShared) {
////    gmp_pool_setup();
//    SchemeParameters params(3, 5, 5, 3, 3, 2);
//    std::default_random_engine rnd(2235);
//    KeysGenerator real(params, rnd);
//    for (int i = 0; i < 10; ++i) {
//      UniformAutomorphismSLPGenerator<int> aut_random(3, &rnd);
//      AutomorphismDescription mock_private_key(3, aut_random);

//      MockParticipant mock(mock_private_key, real.private_key());

//      auto processed_pub_key = mock.process_incoming_public_keys(real.public_keys());

//      auto alice_shared_key = real.make_shared_key(processed_pub_key);
//      auto bob_shared_key = mock.make_shared_key(false);

//      EXPECT_EQ(alice_shared_key().free_reduction(), bob_shared_key().free_reduction());

//      alice_shared_key = mock.make_shared_key(true);
//      bob_shared_key = real.make_shared_key(processed_pub_key, false);

//      EXPECT_EQ(alice_shared_key().free_reduction(), bob_shared_key().free_reduction());
//    }
//  }

//  TEST(FGACryptoTest, SharedKeys) {
////    gmp_pool_setup();
//    SchemeParameters params(3, 5, 4, 2, 2, 1);
//    std::default_random_engine rnd(22345);
////    for (int i = 0; i < 5; ++i) {

////      std::cout << "start" << i << std::endl;
//      KeysGenerator alice(params, rnd);
//      KeysGenerator bob(params, rnd);

//      auto a_pk = alice.public_keys();
//      auto b_pk = bob.public_keys();

//      auto a_processed = bob.process_incoming_public_keys(a_pk);
//      auto b_processed = alice.process_incoming_public_keys(b_pk);

//      auto a_shared_key = alice.make_shared_key(a_processed, true);
//      auto b_shared_key = bob.make_shared_key(b_processed, false);

//      std::cout << "comparing keys" << std::endl;
//      EXPECT_EQ(a_shared_key().free_reduction(), b_shared_key().free_reduction());
////      std::cout << "finish" << i << std::endl;

//      //basic speed for 3, 5, 4, 2, 2, 1 is 6.7s
//      //better 2.9s
////    }
//  }

  TEST(FGACryptoTest, SpeedTest) {
//    gmp_pool_setup();
    SchemeParameters params(3, 5, 4, 25, 25, 50);
    std::default_random_engine rnd(22345);

//    822s for 3, 5, 4, 5, 5, 5 sizes within 3000
//    start line
//    item vert num (val=1548, inv=1524)
//    item vert num (val=1655, inv=1684)
//    item vert num (val=1737, inv=1762)
//    item vert num (val=1855, inv=1929)
//    item vert num (val=1969, inv=1987)
//    line finished
//    key vert num (val=5571, inv=5607)
//    nf key item vert num (val=2099, inv=1962)
//    start line
//    item vert num (val=1418, inv=1483)
//    item vert num (val=1607, inv=1576)
//    item vert num (val=1675, inv=1705)
//    item vert num (val=1861, inv=1904)
//    item vert num (val=1910, inv=1919)
//    line finished
//    key vert num (val=13469, inv=9413)
//    nf key item vert num (val=2215, inv=2116)
//    key vert num (val=13927, inv=9847)
//    nf key item vert num (val=2016, inv=2097)
//    start line
//    item vert num (val=1552, inv=1446)
//    item vert num (val=1630, inv=1587)
//    item vert num (val=1643, inv=1638)
//    item vert num (val=1924, inv=1907)
//    item vert num (val=1907, inv=1836)
//    line finished
//    key vert num (val=13281, inv=9507)
//    nf key item vert num (val=2163, inv=2014)
//    key vert num (val=5957, inv=5844)
//    nf key item vert num (val=2056, inv=2018)
//    almost final key vert num (val=2373, inv=4252)
//    result key item vert num (val=3162, inv=2821)
//    start line
//    item vert num (val=1427, inv=1402)
//    item vert num (val=1643, inv=1610)
//    item vert num (val=1530, inv=1498)
//    item vert num (val=1779, inv=1737)
//    item vert num (val=1766, inv=1698)
//    line finished
//    key vert num (val=4762, inv=4626)
//    nf key item vert num (val=1724, inv=1656)
//    start line
//    item vert num (val=1525, inv=1512)
//    item vert num (val=1684, inv=1641)
//    item vert num (val=1640, inv=1633)
//    item vert num (val=1821, inv=1762)
//    item vert num (val=1727, inv=1674)
//    line finished
//    key vert num (val=11567, inv=7891)
//    nf key item vert num (val=1841, inv=2006)
//    key vert num (val=11936, inv=8618)
//    nf key item vert num (val=1812, inv=1844)
//    key vert num (val=11743, inv=8267)
//    nf key item vert num (val=1898, inv=1907)
//    key vert num (val=5410, inv=5492)
//    nf key item vert num (val=1974, inv=2045)
//    almost final key vert num (val=4212, inv=2418)
//    result key item vert num (val=3162, inv=2821)

    for (int i = 0; i < 1; ++i) {

//      std::cout << "start" << i << std::endl;
      KeysGenerator alice(params, rnd);
      KeysGenerator bob(params, rnd);

      auto a_pk = alice.public_keys();
      auto b_pk = bob.public_keys();

//      std::cout << "s num " << a_pk.s[0].composed_num() << std::endl;
//      std::cout << "private key num " << alice.private_key().composed_num() << std::endl;

      auto a_processed = bob.process_incoming_public_keys(a_pk);
      auto b_processed = alice.process_incoming_public_keys(b_pk);

      auto a_shared_key = alice.make_shared_key(a_processed, true);
      auto b_shared_key = bob.make_shared_key(b_processed, false);

//      std::cout << "a num " << a_shared_key.composed_num() << std::endl;
//      std::cout << "b num " << b_shared_key.composed_num() << std::endl;

//      std::cout << "vertices num=" << slp_vertices_num(a_shared_key()) << std::endl;
//      std::cout << "vertices with unique images num=" << slp_unique_images_length_num(a_shared_key()) << std::endl;
//      std::cout << "height=" << height(a_shared_key()) << std::endl;
//      auto free_red = a_shared_key().free_reduction();
//      std::cout << "free red vertices num=" << slp_vertices_num(free_red) << std::endl;
//      auto normal_form = free_red.normal_form(EndomorphismSLP<int>::NormalFormType::Jez);
//      std::cout << "nf vertices num=" << slp_vertices_num(normal_form) << std::endl;

//      std::cout << "comparing keys" << std::endl;
      EXPECT_EQ(a_shared_key(), b_shared_key());
//      std::cout << "finish" << i << std::endl;

      //basic speed for 3, 5, 4, 2, 2, 1 is 6.7s
      //better 2.9s
      std::cout << std::endl;
    }//3, 5, 4, 5, 5, 5 10s
  }

} // namespace fga_crypto
} // namespace crag
