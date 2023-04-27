#include <gtest/gtest.h>

#include "parallel.h"
#include "walnut_attack.h"

namespace crag {
namespace walnut {
namespace {

using GF256 =
    finitefield::FieldElement<finitefield::IdealGeneratedByPolynomial<finitefield::ZZ<2>, 1, 1, 0, 1, 1, 0, 0, 0, 1>>;

TEST(WalnutAttack, CloakingElementsRemoval1) {
  using FF = GF256;
  using Stabilizer = walnut::StabilizerSquare;
  //  using Stabilizer = walnut::StabilizerDoubleSquare;

  const size_t n = 8;
  const size_t hash_size = 256;
  const size_t seed = 0;

  const auto random_parameters = walnut::randomParameters<FF, Stabilizer>(n, seed).wMaxLength(150);
  const auto random_encoder = walnut::randomEncoder(n, hash_size, seed);
  const auto stabilizer = Stabilizer(random_parameters, 15);
  const auto protocol = walnut::getProtocol(random_parameters, random_encoder, stabilizer);

  const auto private_key = protocol.generatePrivateKey(seed);
  const auto public_key = protocol.computePublicKey(private_key);

  for (size_t i = 0; i < 10; ++i) {
    const auto s1 = protocol.sign(walnut::randomMessageHash(hash_size, seed + i), private_key, seed + i);

    const auto s2 = protocol.sign(walnut::randomMessageHash(hash_size, seed + i + 1), private_key, seed + i + 1);

    const auto without_cloaking_els = removeCloakingElements(protocol, private_key, public_key, s1, s2);

    EXPECT_TRUE(boost::none != without_cloaking_els) << "case " << i << " failed";

    if (without_cloaking_els) {
      EXPECT_TRUE(isUncloakedCorrectly(n, private_key, s1, without_cloaking_els->first)) << "1st case " << i;
      EXPECT_TRUE(isUncloakedCorrectly(n, private_key, s2, without_cloaking_els->second)) << "2nd case " << i;
    }
  }
}

TEST(WalnutAttack, CloakingElementsRemoval2) {
  using FF = GF256;
  using Stabilizer = walnut::StabilizerSquare;
  //  using Stabilizer = walnut::StabilizerDoubleSquare;

  const size_t n = 11;
  const size_t hash_size = 512;
  const size_t seed = 0;

  const auto random_parameters = walnut::randomParameters<FF, Stabilizer>(n, seed).wMaxLength(150);
  const auto random_encoder = walnut::randomEncoder(n, hash_size, seed);
  const auto stabilizer = Stabilizer(random_parameters, 30);
  const auto protocol = walnut::getProtocol(random_parameters, random_encoder, stabilizer);

  const auto private_key = protocol.generatePrivateKey(seed);
  const auto public_key = protocol.computePublicKey(private_key);

  for (size_t i = 0; i < 10; ++i) {
    const auto s1 = protocol.sign(walnut::randomMessageHash(hash_size, seed + i), private_key, seed + i);

    const auto s2 = protocol.sign(walnut::randomMessageHash(hash_size, seed + i + 1), private_key, seed + i + 1);

    const auto without_cloaking_els = removeCloakingElements(protocol, private_key, public_key, s1, s2);

    EXPECT_TRUE(boost::none != without_cloaking_els) << "case " << i << "failed";

    if (without_cloaking_els) {
      EXPECT_TRUE(isUncloakedCorrectly(n, private_key, s1, without_cloaking_els->first)) << "1st case " << i;
      EXPECT_TRUE(isUncloakedCorrectly(n, private_key, s2, without_cloaking_els->second)) << "2nd case " << i;
    }
  }
}

TEST(WalnutAttack, AvailableFlips) {
  std::cout << std::endl << " v1 | v | v2" << std::endl;

  const size_t experiments_count = 100;
  size_t success_count = 0;

  const size_t init_seed = 0;

  for (size_t i = init_seed; i < experiments_count; ++i) {
    const auto protocol = walnut::getProtocolFor256BitsSecurityN11<StabilizerSquare>(i);

    // generate random private key using seed = i
    const auto private_key = protocol.generatePrivateKey(i);
    const auto public_key = protocol.computePublicKey(private_key);

    const auto hash_size = protocol.encoder().hashSize();
    const auto n = protocol.publicParameters().n();
    const auto a = protocol.publicParameters().a();
    const auto b = protocol.publicParameters().b();

    braidgroup::FastIdentityChecker<finitefield::ZZ<10007>> checker(n, i);

    const auto random_msg_hash = randomMessageHash(hash_size, i);
    const auto sig = protocol.sign(random_msg_hash, private_key, i);

    const auto sigma_w1 = public_key.w1Projection().permutation().inverse();

    const auto strand1 = sigma_w1[a];
    const auto strand2 = sigma_w1[b];

    const auto flips = availableFlips<decltype(protocol)::stabilizer_t>(n, strand1, strand2, sig.signature());

    const auto without_v1 = -sig.v1() * sig.signature();
    const auto without_v = sig.v1() * -private_key.w1() * sig.encodedMessageHash() * private_key.w2() * sig.v2();
    const auto without_v2 = sig.signature() * -sig.v2();

    const auto is_without_v1 = parallel::bmap(
        flips, [&](const std::pair<size_t, Word>& f) {
          return !checker.isNonTrivial(-f.second * without_v1);
        });

    const auto is_without_v = parallel::bmap(
        flips, [&](const std::pair<size_t, Word>& f) {
          return !checker.isNonTrivial(-f.second * without_v);
        });

    const auto is_without_v2 = parallel::bmap(
        flips, [&](const std::pair<size_t, Word>& f) {
          return !checker.isNonTrivial(-f.second * without_v2);
        });

    int k = 0;

    if (std::any_of(is_without_v1.begin(), is_without_v1.end(), [](bool b) { return b; })) {
      ++k;
      std::cout << " +  |";
    } else {
      std::cout << " -  |";
    }

    if (std::any_of(is_without_v.begin(), is_without_v.end(), [](bool b) { return b; })) {
      ++k;
      std::cout << " + |";
    } else {
      std::cout << " - |";
    }

    if (std::any_of(is_without_v2.begin(), is_without_v2.end(), [](bool b) { return b; })) {
      ++k;
      std::cout << " +";
    } else {
      std::cout << " -";
    }

    if (k == 3) {
      ++success_count;
    }

    std::cout << std::endl;
  }

  std::cout << "Success: " << success_count << "/" << experiments_count << std::endl;
}

TEST(WalnutAttack, Ex_128_bits_security) {
  const size_t e = 0; // init seed

  const auto protocol = walnut::getProtocolFor128BitsSecurity(e);

  // generate random private key using seed = e
  const auto private_key = protocol.generatePrivateKey(e);
  const auto public_key = protocol.computePublicKey(private_key);

  std::cout << std::endl;

  std::cout << "|w1|  = " << private_key.w1().length() << ", ";
  std::cout << "|w2|  = " << private_key.w2().length() << std::endl;

  bool success = false;
  const size_t attempt_number = 10;

  for (size_t attempt = 0; (attempt < attempt_number) && !success; ++attempt) {
    std::cout << "Attempt #" << attempt << std::endl;

    std::mt19937_64 g(attempt);

    const auto fake_private_key = walnut::attack(protocol, private_key, public_key, g);

    if (success = (fake_private_key != boost::none)) {
      if (success &= walnut::checkFakePrivateKey(protocol, *fake_private_key, private_key, public_key, g)) {
        std::cout << "Computed private key is GOOD." << std::endl;
      } else {
        std::cout << "Computed private key is BAD." << std::endl;
      }
    }
  }
}

TEST(WalnutAttack, Ex_256_bits_security) {
  const size_t e = 0; // init seed

  const auto protocol = walnut::getProtocolFor256BitsSecurity(e);

  // generate random private key using seed = e
  const auto private_key = protocol.generatePrivateKey(e);
  const auto public_key = protocol.computePublicKey(private_key);

  std::cout << std::endl;

  std::cout << "|w1|  = " << private_key.w1().length() << ", ";
  std::cout << "|w2|  = " << private_key.w2().length() << std::endl;

  bool success = false;
  const size_t attempt_number = 10;

  for (size_t attempt = 0; (attempt < attempt_number) && !success; ++attempt) {
    std::cout << "Attempt #" << attempt << std::endl;

    std::mt19937_64 g(attempt);

    const auto fake_private_key = walnut::attack(protocol, private_key, public_key, g);

    if (success = (fake_private_key != boost::none)) {
      if (success &= walnut::checkFakePrivateKey(protocol, *fake_private_key, private_key, public_key, g)) {
        std::cout << "Computed private key is GOOD." << std::endl;
      } else {
        std::cout << "Computed private key is BAD." << std::endl;
      }
    }
  }
}

TEST(WalnutAttack, Ex_256_bits_security_n_11) {
  const size_t e = 0; // init seed

  const auto protocol = walnut::getProtocolFor256BitsSecurityN11(e);

  // generate random private key using seed = e
  const auto private_key = protocol.generatePrivateKey(e);
  const auto public_key = protocol.computePublicKey(private_key);

  std::cout << std::endl;

  std::cout << "|w1|  = " << private_key.w1().length() << ", ";
  std::cout << "|w2|  = " << private_key.w2().length() << std::endl;

  bool success = false;
  const size_t attempt_number = 10;

  for (size_t attempt = 0; (attempt < attempt_number) && !success; ++attempt) {
    std::cout << "Attempt #" << attempt << std::endl;

    std::mt19937_64 g(attempt);

    const auto fake_private_key = walnut::attack(protocol, private_key, public_key, g);

    if (success = (fake_private_key != boost::none)) {
      if (success &= walnut::checkFakePrivateKey(protocol, *fake_private_key, private_key, public_key, g)) {
        std::cout << "Computed private key is GOOD." << std::endl;
      } else {
        std::cout << "Computed private key is BAD." << std::endl;
      }
    }
  }
}
} // namespace
} // namespace walnut
} // namespace crag
