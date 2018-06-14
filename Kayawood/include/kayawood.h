#pragma once

#ifndef CRAG_KAYAWOOD_H
#define CRAG_KAYAWOOD_H

#include "ShortBraidForm.h"
#include "Word.h"
#include "cloaking_element.h"
#include "colored_burau.h"
#include "kayawood_cloaking.h"
#include "kayawood_parameters.h"
#include "parallel.h"
#include "random_word.h"
#include "stochastic_rewrite.h"

namespace crag {
namespace kayawood {

//! Default rewriting algorithm used for braid obfuscation in Kayawood.
struct DehornoyObfuscator {
  Word operator()(size_t n, const Word& w) const {
    return shortenBraid2(n, w);
  }
};

//! First applies Right Garside normal form and then Dehornoy reduction.
struct GarsideDehornoyObfuscator {
  Word operator()(size_t n, const Word& w) const {
    return shortBraidForm(n, w);
  }
};

//! For testing purposes.
struct TrivialObfuscator {
  Word operator()(size_t n, const Word& w) const {
    return w;
  }
};

//! Obfuscator applying stochastic rewrite and then Dehornoy reduction.
class StochasticRewriteObfuscator {
public:
  StochasticRewriteObfuscator(
      std::vector<size_t> partition, size_t min_block_size, size_t max_block_size, size_t iter_num, size_t seed)
      : partition_(std::move(partition))
      , min_block_size_(min_block_size)
      , max_block_size_(max_block_size)
      , iter_num_(iter_num)
      , seed_(seed) {}

  Word operator()(size_t n, const Word& w) const {
    static std::mt19937_64 g(seed_);
    const auto stochastic_w =
        stochasticrewrite::stochasticRewrite(w, partition_, min_block_size_, max_block_size_, iter_num_, g);
    return stochastic_w;
  }

private:
  std::vector<size_t> partition_;

  size_t min_block_size_;
  size_t max_block_size_;
  size_t iter_num_;

  size_t seed_;
};

StochasticRewriteObfuscator getStochasticRewriteObfuscator(size_t n, size_t seed);

using coloredburau::CBProjectionElement;

template <typename T>
class ProtocolInstance {
public:
  // clang-format off
  ProtocolInstance(
      PublicParameters<T> parameters,
      std::vector<Word> betas,
      std::vector<Word> betas_conjugates,
      Word z,
      Word alpha,
      Word alice_private_key,
      Word bob_private_key,
      Word alice_public_key,
      Word bob_public_key,
      CBProjectionElement<T> shared_key)
      : parameters_(std::move(parameters))
      , betas_(std::move(betas))
      , betas_conjugates_(std::move(betas_conjugates))
      , z_(std::move(z))
      , alpha_(std::move(alpha))
      , alice_private_key_(std::move(alice_private_key))
      , bob_private_key_(std::move(bob_private_key))
      , alice_public_key_(std::move(alice_public_key))
      , bob_public_key_(std::move(bob_public_key))
      , shared_key_(std::move(shared_key)) {}
  // clang-format on

  const PublicParameters<T>& parameters() const {
    return parameters_;
  }

  const std::vector<Word>& betas() const {
    return betas_;
  }

  const std::vector<Word>& betasConjugates() const {
    return betas_conjugates_;
  }

  const Word& z() const {
    return z_;
  }

  const Word& alpha() const {
    return alpha_;
  }

  const Word& alicePrivateKey() const {
    return alice_private_key_;
  }

  const Word& bobPrivateKey() const {
    return bob_private_key_;
  }

  const Word& alicePublicKey() const {
    return alice_public_key_;
  }

  const Word& bobPublicKey() const {
    return bob_public_key_;
  }

  const CBProjectionElement<T>& sharedKey() const {
    return shared_key_;
  }

private:
  PublicParameters<T> parameters_;

  std::vector<Word> betas_;
  std::vector<Word> betas_conjugates_;

  Word z_;
  Word alpha_;

  Word alice_private_key_;
  Word bob_private_key_;

  Word alice_public_key_;
  Word bob_public_key_;

  CBProjectionElement<T> shared_key_;
};

//! Describes the Kayawood protocol.
//! Parameters:
//! 1. Finite field F containing at least 32 elements, specified by T.
//! 2. Rewriting algorithms for braid words - callable object with Word operator()(size_t n, const Word& w).
template <
    typename T,
    typename BetasObfuscator = GarsideDehornoyObfuscator,
    typename PublicKeyObfuscator = DehornoyObfuscator,
    typename Stabilizer = StabilizerSquare>
class Protocol {
public:
  using field_t = T;
  using stabilizer_t = Stabilizer;

  Protocol(PublicParameters<T> parameters, Stabilizer stabilizer)
      : parameters_(std::move(parameters))
      , stabilizer_(std::move(stabilizer)) {}

  Protocol(
      PublicParameters<T> parameters,
      BetasObfuscator betas_obfuscator,
      PublicKeyObfuscator public_key_obfuscator,
      Stabilizer stabilizer)
      : parameters_(std::move(parameters))
      , betas_obfuscator_(betas_obfuscator)
      , public_key_obfuscator_(public_key_obfuscator)
      , stabilizer_(std::move(stabilizer)) {}

  template <typename URNG>
  ProtocolInstance<T> generateInstance(URNG& g) const {
    const auto n = parameters_.n();

    auto z = generateZ_(n, parameters_.zMinLength(), parameters_.zMaxLength(), g);

    auto betas = generateBetas_(n, parameters_.r(), parameters_.betaMinLength(), parameters_.betaMaxLength(), g);
    //    std::cout << "Started beta conjugates" << std::endl;
    auto betas_conjugates = computeBetasConjugates_(n, betas, z);
    //    std::cout << "Finished beta conjugates" << std::endl;

    auto alpha = generateAlpha_(n, parameters_.alphaMinLength(), parameters_.alphaMaxLength(), g);

    auto alice_private_key = z * alpha * -z;
    auto sigma_a = coloredburau::permutation(n, alice_private_key);

    auto bob_private_key = generateBobPrivateKey_(
        n, parameters_.bobPrivateKeyMinLength(), parameters_.bobPrivateKeyMaxLength(), betas_conjugates, g);

    auto sigma_b = coloredburau::permutation(n, bob_private_key);

    const auto pub_b = stabilizer_(sigma_a, bob_private_key, g);
    //    std::cout << "Started obfuscation of pub_B, |pub_B| = " << pub_b.size() << std::endl;
    auto bob_public_key = public_key_obfuscator_(n, pub_b);
    //    std::cout << "Finished obfuscation of pub_B, |pub_B| = " << bob_public_key.size() << std::endl;

    const auto pub_a = stabilizer_(sigma_b, alice_private_key, g);
    //    std::cout << "Started obfuscation of pub_A, |pub_A| = " << pub_a.size() << std::endl;
    auto alice_public_key = public_key_obfuscator_(n, pub_a);
    //    std::cout << "Finished obfuscation of pub_A, |pub_A| = " << alice_public_key.size() << std::endl;

    const auto alice_shared_key = coloredburau::project(alice_private_key * bob_public_key, parameters_.tValues());
    const auto bob_shared_key = coloredburau::project(bob_private_key * alice_public_key, parameters_.tValues());

    // true_shared_key is the expected value of the shared key
    auto true_shared_key = coloredburau::project(alice_private_key * bob_private_key, parameters_.tValues());

    if ((bob_shared_key != true_shared_key) || (alice_shared_key != true_shared_key)) {
      throw std::logic_error("Shared keys are not equal.");
    }

    // clang-format off
    return ProtocolInstance<T>(
        parameters_,
        std::move(betas),
        std::move(betas_conjugates),
        std::move(z),
        std::move(alpha),
        std::move(alice_private_key),
        std::move(bob_private_key),
        std::move(alice_public_key),
        std::move(bob_public_key),
        std::move(true_shared_key)
    );
    // clang-format on
  }

  ProtocolInstance<T> generateInstance(size_t seed) const {
    std::mt19937_64 g(seed);
    return generateInstance(g);
  }

  const PublicParameters<T>& parameters() const {
    return parameters_;
  }

private:
  PublicParameters<T> parameters_;

  BetasObfuscator betas_obfuscator_;
  PublicKeyObfuscator public_key_obfuscator_;

  Stabilizer stabilizer_;

  template <typename URNG>
  Word generateAlpha_(size_t n, size_t min_length, size_t max_length, URNG& g) const {
    return random::randomWord(n / 2 - 1, min_length, max_length, g);
  }

  //! Generates beta_1,...,beta_r from the subgroup generated by b_{n/2 + 1},...,b_{n-1} in B_n
  template <typename URNG>
  std::vector<Word> generateBetas_(size_t n, size_t r, size_t min_length, size_t max_length, URNG& g) const {
    const auto generators_count = n / 2 - 1;

    std::vector<Word> bs;
    bs.reserve(generators_count);

    for (size_t i = 0; i < generators_count; ++i) {
      bs.emplace_back(generators_count + 2 + i);
    }

    std::vector<Word> betas;
    betas.reserve(r);

    for (size_t i = 0; i < r; ++i) {
      const auto random_word = random::randomWord(generators_count, min_length, max_length, g);
      betas.push_back(random_word.replaceGenerators(bs));
    }

    return betas;
  }

  std::vector<Word> computeBetasConjugates_(size_t n, const std::vector<Word>& betas, const Word& z) const {
    return parallel::map(betas, [&, this](const Word& beta) { return betas_obfuscator_(n, z * beta * -z); });
  }

  template <typename URNG>
  Word generateZ_(size_t n, size_t min_length, size_t max_length, URNG& g) const {
    const size_t max_attempts = 100;

    const size_t mixed_strands_goal = n / 4;

    for (size_t i = 0; i < max_attempts; ++i) {
      const auto z = random::randomWord(n - 1, min_length, max_length, g);

      // Check that z satisfies the z-condition: mixed_strands = mixed_strands_goal
      const auto p = coloredburau::permutation(n, z);

      size_t mixed_strands = 0;

      for (size_t j = 0; j < n / 2; ++j) {
        if (p[j] >= n / 2) {
          ++mixed_strands;
        }
      }

      if (mixed_strands == mixed_strands_goal) {
        return z;
      }
    }

    throw std::runtime_error("Cannot generate a proper z braid.");
  }

  template <typename URNG>
  Word generateBobPrivateKey_(
      size_t n, size_t min_length, size_t max_length, const std::vector<Word>& betas_conjugates, URNG& g) const {
    const auto r = betas_conjugates.size();

    return random::randomWord(r, min_length, max_length, g).replaceGenerators(betas_conjugates);
  }
};

//! Returns true iff shared key can be obtained from public key.
template <typename T>
bool isBadInstance(const ProtocolInstance<T>& instance) {
  const auto& t_values = instance.parameters().tValues();
  const auto& pub_a = instance.alicePublicKey();
  const auto& pub_b = instance.bobPublicKey();

  const auto shared_key_from_pub_a_pub_b = coloredburau::project(pub_a * pub_b, t_values);

  if (instance.sharedKey() == shared_key_from_pub_a_pub_b) {
    return true;
  }

  const auto shared_key_from_pub_b_pub_a = coloredburau::project(pub_b * pub_a, t_values);

  if (instance.sharedKey() == shared_key_from_pub_b_pub_a) {
    return true;
  }

  return false;
}

template <typename T, typename BetasObfuscator, typename PublicKeyObfuscator, typename Stabilizer = StabilizerSquare>
Protocol<T, BetasObfuscator, PublicKeyObfuscator, Stabilizer> getProtocol(
    PublicParameters<T> parameters,
    BetasObfuscator betas_obfuscator,
    PublicKeyObfuscator public_key_obfuscator,
    Stabilizer stabilizer) {
  return Protocol<T, BetasObfuscator, PublicKeyObfuscator, Stabilizer>(
      std::move(parameters), std::move(betas_obfuscator), std::move(public_key_obfuscator), std::move(stabilizer));
}

template <
    typename T,
    typename BetasObfuscator = GarsideDehornoyObfuscator,
    typename PublicKeyObfuscator = DehornoyObfuscator,
    typename Stabilizer = StabilizerSquare>
Protocol<T, BetasObfuscator, PublicKeyObfuscator, Stabilizer>
getProtocol(PublicParameters<T> parameters, Stabilizer stabilizer) {
  return getProtocol(std::move(parameters), BetasObfuscator(), PublicKeyObfuscator(), std::move(stabilizer));
}
} // namespace kayawood
} // namespace crag

#endif // CRAG_KAYAWOOD_H
