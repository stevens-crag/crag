#pragma once

#ifndef CRAG_WALNUT_H
#define CRAG_WALNUT_H

#include <vector>

#include "ShortBraidForm.h"
#include "Word.h"
#include "colored_burau.h"
#include "random_word.h"
#include "walnut_cloaking.h"
#include "walnut_encoding.h"
#include "walnut_parameters.h"

namespace crag {
namespace walnut {

using coloredburau::CBProjectionElement;

//! Default rewriting algorithm used for signature obfuscation in Walnut.
//! First applies Right Garside normal form and then Dehornoy reduction.
struct GarsideDehornoyObfuscator {
  Word operator()(size_t n, const Word& w) const {
    return shortBraidForm(n, w);
  }
};

//! Private key for Walnut - a pair of freely reduced braids (w1, w2)
//! such that w1, w2, and w2*w1 are not in the pure braid subgroup.
class PrivateKey {
public:
  PrivateKey(Word w1, Word w2)
      : w1_(std::move(w1))
      , w2_(std::move(w2)) {}

  const Word& w1() const {
    return w1_;
  }

  const Word& w2() const {
    return w2_;
  }

private:
  Word w1_;
  Word w2_;
};

//! Public key for Walnut - projection of a private key.
template <typename T>
class PublicKey {
public:
  PublicKey(CBProjectionElement<T> w1_projection, CBProjectionElement<T> w2_projection)
      : w1_projection_(std::move(w1_projection))
      , w2_projection_(std::move(w2_projection)) {
    if (w1_projection_.tValues() != w2_projection_.tValues()) {
      throw std::invalid_argument("Elements must use the same t-values.");
    }
  }

  const CBProjectionElement<T>& w1Projection() const {
    return w1_projection_;
  }

  const CBProjectionElement<T>& w2Projection() const {
    return w2_projection_;
  }

private:
  CBProjectionElement<T> w1_projection_;
  CBProjectionElement<T> w2_projection_;
};

//! Signature of a message in Walnut.
class Signature {
public:
  Signature(msg_hash_t msg_hash, Word encoded_msg_hash, Word signature, Word v, Word v1, Word v2)
      : msg_hash_(std::move(msg_hash))
      , encoded_msg_hash_(std::move(encoded_msg_hash))
      , signature_(std::move(signature))
      , v_(std::move(v))
      , v1_(std::move(v1))
      , v2_(std::move(v2)) {}

  const msg_hash_t& messageHash() const {
    return msg_hash_;
  }

  const Word& encodedMessageHash() const {
    return encoded_msg_hash_;
  }

  const Word& signature() const {
    return signature_;
  }

  const Word& v() const {
    return v_;
  }

  const Word& v1() const {
    return v1_;
  }

  const Word& v2() const {
    return v2_;
  }

private:
  msg_hash_t msg_hash_;
  Word encoded_msg_hash_;
  Word signature_;

  Word v_;
  Word v1_;
  Word v2_;
};

template <
    typename T,
    typename Obfuscator = GarsideDehornoyObfuscator,
    typename Encoder = DefaultEncoder,
    typename Stabilizer = StabilizerSquare>
class Protocol {
public:
  using obfuscator_t = Obfuscator;
  using encoder_t = Encoder;
  using stabilizer_t = Stabilizer;

  Protocol(PublicParameters<T> parameters, Encoder encoder, Stabilizer stabilizer)
      : parameters_(std::move(parameters))
      , encoder_(std::move(encoder))
      , stabilizer_(std::move(stabilizer)) {}

  template <typename URNG>
  PrivateKey generatePrivateKey(URNG& g) const {
    const auto n = parameters_.n();
    const auto min_length = parameters_.wMinLength();
    const auto max_length = parameters_.wMaxLength();

    Word w1, w2;

    // Generate w1, w2 such that w1, w2, and w2*w1 are not in the pure braid subgroup.
    // With high probability this is the case, so use while(true) loop.
    while (true) {
      w1 = random::randomWord(n - 1, min_length, max_length, g);
      w2 = random::randomWord(n - 1, min_length, max_length, g);

      const auto sigma_w1 = coloredburau::permutation(n, w1);

      if (sigma_w1.isTrivial()) {
        continue;
      }

      const auto sigma_w2 = coloredburau::permutation(n, w1);

      if (sigma_w2.isTrivial()) {
        continue;
      }

      if ((sigma_w1 * sigma_w2).isTrivial()) {
        continue;
      }

      return PrivateKey(std::move(w1), std::move(w2));
    }
  }

  PrivateKey generatePrivateKey(size_t seed) const {
    std::mt19937_64 g(seed);
    return generatePrivateKey(g);
  }

  PublicKey<T> computePublicKey(const PrivateKey& private_key) const {
    const auto& t_values = parameters_.tValues();

    return PublicKey<T>(
        coloredburau::project(private_key.w1(), t_values), coloredburau::project(private_key.w2(), t_values));
  }

  template <typename URNG>
  Signature sign(const msg_hash_t& msg_hash, const PrivateKey& private_key, URNG& g) const {
    const auto n = parameters_.n();

    const auto& w1 = private_key.w1();
    const auto& w2 = private_key.w2();

    const auto sigma_w1 = coloredburau::permutation(n, w1);
    const auto sigma_w2 = coloredburau::permutation(n, w2);

    auto v = stabilizer_(Permutation(n), g);
    auto v1 = stabilizer_(sigma_w1, g);
    auto v2 = stabilizer_(sigma_w2, g);

    auto encoded_hash = encoder_(msg_hash);

    const auto signature = v1 * -w1 * v * encoded_hash * w2 * v2;

    return Signature(
        msg_hash, std::move(encoded_hash), obfuscator_(n, signature), std::move(v), std::move(v1), std::move(v2));
  }

  Signature sign(const msg_hash_t& msg_hash, const PrivateKey& private_key, size_t seed) const {
    std::mt19937_64 g(seed);
    return sign(msg_hash, private_key, g);
  }

  Signature signWithoutCloakingEls(const msg_hash_t& msg_hash, const PrivateKey& private_key) const {
    auto encoded_hash = encoder_(msg_hash);

    const auto& w1 = private_key.w1();
    const auto& w2 = private_key.w2();

    const auto signature = -w1 * encoded_hash * w2;

    return Signature(msg_hash, std::move(encoded_hash), signature, Word(), Word(), Word());
  }

  bool verify(const Signature& signature, const PublicKey<T>& public_key) const {
    return verify(signature.messageHash(), signature.signature(), public_key);
  }

  bool verify(const msg_hash_t& msg_hash, const Word& signature, const PublicKey<T>& public_key) const {
    if (public_key.w1Projection().tValues() != parameters_.tValues()) {
      throw std::invalid_argument("T-values of the public key don't match t-values of the protocol.");
    }

    const auto encoded_hash = encoder_(msg_hash);

    const auto lhs = (public_key.w1Projection() * signature).matrix();

    const auto rhs =
        coloredburau::project(encoded_hash, parameters_.tValues()).matrix() * public_key.w2Projection().matrix();

    return lhs == rhs;
  }

  const PublicParameters<T>& publicParameters() const {
    return parameters_;
  }

  const Encoder& encoder() const {
    return encoder_;
  }

  const Stabilizer& stabilizer() const {
    return stabilizer_;
  }

private:
  PublicParameters<T> parameters_;
  Encoder encoder_;
  Stabilizer stabilizer_;

  Obfuscator obfuscator_;
};

template <
    typename T,
    typename Obfuscator = GarsideDehornoyObfuscator,
    typename Encoder = DefaultEncoder,
    typename Stabilizer = StabilizerSquare>
Protocol<T, Obfuscator, Encoder, Stabilizer>
getProtocol(PublicParameters<T> parameters, Encoder encoder, Stabilizer stabilizer) {
  return Protocol<T, Obfuscator, Encoder, Stabilizer>(std::move(parameters), std::move(encoder), std::move(stabilizer));
};
} // namespace walnut
} // namespace crag

#endif // CRAG_WALNUT_H
