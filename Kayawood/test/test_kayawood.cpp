#include <gtest/gtest.h>

#include "kayawood.h"
#include "random_permutation.h"

namespace crag {
namespace kayawood {
namespace {

using ZZ5 = finitefield::ZZ<5>;
using GF32 = finitefield::FieldElement<finitefield::IdealGeneratedByPolynomial<finitefield::ZZ<2>, 1, 0, 1, 0, 0, 1>>;
using GF256 =
    finitefield::FieldElement<finitefield::IdealGeneratedByPolynomial<finitefield::ZZ<2>, 1, 1, 0, 1, 1, 0, 0, 0, 1>>;

TEST(Kayawood, Ex_01) {
  // clang-format off
  const auto parameters = PublicParameters<ZZ5>(16, 2, 4, {
      ZZ5(3), ZZ5(2), ZZ5(1), ZZ5(4), ZZ5(1), ZZ5(2), ZZ5(3), ZZ5(2),
      ZZ5(3), ZZ5(2), ZZ5(4), ZZ5(4), ZZ5(2), ZZ5(2), ZZ5(3), ZZ5(3),
  });
  // clang-format on
  const auto stabilizer = StabilizerSquare(parameters, 15);

  const auto protocol = Protocol<ZZ5>(parameters, stabilizer);

  std::mt19937 g(0);

  const auto instance = protocol.generateInstance(g);
}

TEST(Kayawood, RandomParams1) {
  using FF = GF32;
  const size_t n = 16;
  const size_t L = 15;

  const auto random_params = randomParameters<FF, StabilizerSquare>(n, 0);
  const auto stabilizer = StabilizerSquare(random_params, L);

  const auto protocol = getProtocol(random_params, stabilizer);

  const auto instance = protocol.generateInstance(0);
}

TEST(Kayawood, RandomParams128bit) {
  using FF = GF32;
  const size_t n = 16;
  const size_t L = 15;
  const size_t r = 32;
  const size_t bob_private_key_length = 22;

  std::cout << std::endl;

  size_t good_instances = 0;

  for (size_t i = 0; i < 100; ++i) {
    const auto random_parameters = randomParameters<FF, StabilizerSquare>(n, i)
                                       .r(r)
                                       .zMinLength(180)
                                       .zMaxLength(250)
                                       .alphaMinLength(300)
                                       .alphaMaxLength(400)
                                       .betaMinLength(50)
                                       .betaMaxLength(100)
                                       .bobPrivateKeyMinLength(bob_private_key_length)
                                       .bobPrivateKeyMaxLength(bob_private_key_length);

    const auto stabilizer = StabilizerSquare(random_parameters, L);
    const auto protocol = getProtocol(random_parameters, TrivialObfuscator(), TrivialObfuscator(), stabilizer);

    const auto instance = protocol.generateInstance(i);

    const auto is_bad_instance = isBadInstance(instance);

    if (!is_bad_instance) {
      ++good_instances;
    }

    std::cout << "==============================================================================" << std::endl;
    std::cout << "Bad hared key for instance #" << i << ": " << std::boolalpha << is_bad_instance << std::endl;

    const auto n = protocol.parameters().n();
    const auto a = protocol.parameters().a();
    const auto b = protocol.parameters().b();

    const auto sigma_a = coloredburau::permutation(n, instance.alicePrivateKey());
    const auto sigma_b = coloredburau::permutation(n, instance.bobPrivateKey());
    const auto sigma_z = coloredburau::permutation(n, instance.z());

    if (is_bad_instance) {
      EXPECT_TRUE((sigma_a[a] == a && sigma_a[b] == b) || (sigma_b[a] == a && sigma_b[b] == b));
    }

    std::cout << "a = " << a << ", b = " << b << std::endl;
    std::cout << "sigma_a = " << sigma_a << std::endl;
    std::cout << "sigma_b = " << sigma_b << std::endl;
    std::cout << "sigma_z = " << sigma_z << std::endl;

    std::cout << "sigma_alpha = " << (sigma_z.inverse() * sigma_a * sigma_z) << std::endl;
    std::cout << "sigma_beta = " << (sigma_z.inverse() * sigma_b * sigma_z) << std::endl;

    const auto bad_z_images = (sigma_z[a] < n / 2) == (sigma_z[b] < n / 2);
    std::cout << "bad z images: " << std::boolalpha << bad_z_images << std::endl;

    EXPECT_TRUE(!bad_z_images || bad_z_images && is_bad_instance);
  }

  std::cout << "Good instances: " << good_instances << std::endl;
}

TEST(Kayawood, RandomParams256bit) {
  using FF = GF256;
  const size_t n = 16;
  const size_t L = 30;
  const size_t r = 32;
  const size_t bob_private_key_length = 43;

  std::cout << std::endl;

  size_t good_instances = 0;

  for (size_t i = 0; i < 100; ++i) {
    const auto random_parameters = randomParameters<FF, StabilizerSquare>(n, i)
                                       .r(r)
                                       .zMinLength(300)
                                       .zMaxLength(400)
                                       .alphaMinLength(300)
                                       .alphaMaxLength(400)
                                       .betaMinLength(100)
                                       .betaMaxLength(200)
                                       .bobPrivateKeyMinLength(bob_private_key_length)
                                       .bobPrivateKeyMaxLength(bob_private_key_length);

    const auto stabilizer = StabilizerSquare(random_parameters, L);
    const auto protocol = getProtocol(random_parameters, TrivialObfuscator(), TrivialObfuscator(), stabilizer);

    const auto instance = protocol.generateInstance(i);

    //    std::cout << "|priv_A| = " << shortenBraid(n, instance.alicePrivateKey()).size() << std::endl;
    //
    //    std::cout << "|beta_i^z| = ";
    //    for (const auto& w : instance.betasConjugates()) {
    //      std::cout << w.size() << ", ";
    //    }
    //    std::cout << std::endl;

    const auto is_bad_instance = isBadInstance(instance);

    if (!is_bad_instance) {
      ++good_instances;
    }

    std::cout << "Bad hared key for instance #" << i << ": " << std::boolalpha << is_bad_instance << std::endl;
  }

  std::cout << "Good instances: " << good_instances << std::endl;
}

TEST(Kayawood, RandomParams256bitMultipleCloaking) {
  using FF = GF256;
  const size_t n = 16;
  const size_t r = 32;
  const size_t bob_private_key_length = 43;

  std::cout << std::endl;

  size_t good_instances = 0;

  for (size_t i = 0; i < 100; ++i) {
    const auto random_parameters = randomParameters<FF, StabilizerManyShort>(n, i)
                                       .r(r)
                                       .zMinLength(300)
                                       .zMaxLength(400)
                                       .alphaMinLength(300)
                                       .alphaMaxLength(400)
                                       .betaMinLength(100)
                                       .betaMaxLength(200)
                                       .bobPrivateKeyMinLength(bob_private_key_length)
                                       .bobPrivateKeyMaxLength(bob_private_key_length);

    const auto stabilizer = StabilizerManyShort(random_parameters, 30, 50, 30);
    const auto protocol = getProtocol(random_parameters, TrivialObfuscator(), TrivialObfuscator(), stabilizer);

    const auto instance = protocol.generateInstance(i);

    const auto is_bad_instance = isBadInstance(instance);

    if (!is_bad_instance) {
      ++good_instances;
    }

    std::cout << "Bad hared key for instance #" << i << ": " << std::boolalpha << is_bad_instance << std::endl;
  }

  std::cout << "Good instances: " << good_instances << std::endl;
}

TEST(Kayawood, BadSharedKeyWithoutCloakingElements) {
  using FF = GF256;
  const auto n = 16;

  const auto parameters = randomParameters<FF, StabilizerTrivial>(n, 0);
  const auto protocol = getProtocol(parameters, StabilizerTrivial());

  const auto instance = protocol.generateInstance(0);

  EXPECT_TRUE(isBadInstance(instance));
}

TEST(Kayawood, StochasticRewrite) {
  using FF = GF256;

  const size_t n = 16;
  const size_t L = 30;

  const auto parameters = randomParameters<FF, StabilizerSquare>(n, 0);
  const auto stabilizer = StabilizerSquare(parameters, L);

  const auto protocol =
      getProtocol(parameters, GarsideDehornoyObfuscator(), getStochasticRewriteObfuscator(n, 1), stabilizer);

  const auto instance = protocol.generateInstance(0);
}

TEST(Kayawood, Stabilizers) {
  using FF = GF256;

  const size_t n = 16;

  const auto parameters = randomParameters<FF, StabilizerManyShort>(n, 0);
  const auto stabilizer_square = StabilizerSquare(parameters, 0);
  const auto stabilizer_shortmany = StabilizerManyShort(parameters, 30, 50, 30);

  std::mt19937_64 g(0);

  for (size_t i = 0; i < 20; ++i) {
    const auto sigma = random::randomPermutation(n, g);
    //    const auto cloaking_el = stabilizer_square(sigma, Word(), g);
    const auto cloaking_el = stabilizer_shortmany(sigma, Word(), g);

    std::cout << cloaking_el.size() << std::endl;
  }
}
} // namespace
} // namespace kayawood
} // namespace crag
