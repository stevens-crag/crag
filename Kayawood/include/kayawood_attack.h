#pragma once

#ifndef CRAG_KAYAWOOD_ATTACK_H
#define CRAG_KAYAWOOD_ATTACK_H

#include <fstream>
#include <future>

#include "kayawood.h"

#include "LinkedBraidStructure.h"
#include "ThLeftNormalForm.h"
#include "braid_group.h"
#include "fast_identity_check.h"

namespace crag {
namespace kayawood {

//! Checks if u \in B_n commutes with each w in words.
static bool doesCommuteWithTuple(size_t n, const Word& u, const std::vector<Word>& words) {
  static const crag::braidgroup::FastIdentityChecker<finitefield::ZZ<1237>> checker(n);

  for (const auto& w : words) {
    const auto commutator = w * u * -w * -u;

    // quick fast check first
    if (checker.isNonTrivial(commutator)) {
      return false;
    }

    // complete check
    if (!isTrivialBraid(n, commutator)) {
      return false;
    }
  }

  return true;
}

//! Computes the total length of all words.
static size_t totalLength(const std::vector<Word>& words) {
  size_t result = 0;

  for (const auto& w : words) {
    result += w.length();
  }

  return result;
}

//! For a vector of words returns a vector of shortened conjugates.
static std::vector<Word> conjugate(size_t n, const std::vector<Word>& words, const Word& c) {
  std::vector<Word> conjugated_words;
  conjugated_words.reserve(words.size());

  for (const auto& w : words) {
    conjugated_words.push_back(shortenBraid(n, (-c) * w * c));
  }

  return conjugated_words;
}

//! Tries to recover (z, betas) from betas^{z^{-1}}.
static std::pair<Word, std::vector<Word>> recoverZ(size_t n, const std::vector<Word>& betas_conjugates) {
  // Prepare the conjugators
  std::vector<Word> conjugators;
  conjugators.reserve(2 * n);

  for (int i = 1; i < n; ++i) {
    conjugators.push_back(Word(i));
    conjugators.push_back(Word(-i));
  }

  const auto initial_words = std::make_pair(totalLength(betas_conjugates), betas_conjugates);

  std::map<std::pair<size_t, std::vector<Word>>, Word> checked_elements;

  std::map<std::pair<size_t, std::vector<Word>>, Word> unchecked_elements = {
      {initial_words, Word()},
  };

  size_t best_len = initial_words.first;

  std::cout << "Beta total length: ";

  for (size_t attempt = 0; attempt < 20 && !unchecked_elements.empty(); ++attempt) {
    const auto current = *unchecked_elements.begin();

    unchecked_elements.erase(current.first);
    checked_elements.emplace(current);

    const auto& cur_tuple = current.first.second;
    const auto& cur_z = current.second;
    const auto cur_len = current.first.first;

    std::cout << cur_len << ",";

    // add additional conjugator
    conjugators.push_back(cur_tuple[0].subword(0, 50));

    // map each conjugator c to cur_tuple^c
    const auto conjugated_tuples =
        parallel::map<Word, std::vector<Word>>(conjugators, [&](const Word& c) { return conjugate(n, cur_tuple, c); });

    for (size_t i = 0; i < conjugators.size(); ++i) {
      const auto& c = conjugators[i];
      const auto& new_tuple = conjugated_tuples[i];
      const auto new_total_length = totalLength(new_tuple);

      const auto pr = std::make_pair(new_total_length, new_tuple);

      if ((unchecked_elements.count(pr) == 0) && (checked_elements.count(pr) == 0)) {
        unchecked_elements.emplace(pr, cur_z * c);

        if (new_total_length < best_len) {
          best_len = new_total_length;
          attempt = 0;
        }
      }
    }

    // remove additional conjugator
    conjugators.pop_back();
  }

  std::cout << std::endl;

  const auto best = *checked_elements.begin();

  return std::make_pair(best.second, best.first.second);
}


static Word dropPairs(size_t n, size_t a, size_t b, const Word& given) {
  Word w = given;

  if (a > b) {
    std::swap(a, b);
  }

  size_t total_drop = 0;

  for (bool progress = true; progress;) {
    progress = false;

    Permutation perm1(n), perm2(n);

    auto w_as_vec = w.toVector();

    for (size_t pos = 0; pos + 1 < w_as_vec.size(); ++pos) {
      // Compute strands
      const auto ind = std::abs(w_as_vec[pos]);

      auto strand1 = perm2[ind - 1];
      auto strand2 = perm2[ind];

      perm1.change(strand1, strand2);
      perm2.change(ind - 1, ind);

      if (strand1 > strand2) {
        std::swap(strand1, strand2);
      }

      // Drop pair
      if (strand1 == a && strand2 == b && w_as_vec[pos] == w_as_vec[pos + 1]) {
        w_as_vec.erase(w_as_vec.begin() + pos, w_as_vec.begin() + pos + 2);
        total_drop++;
        progress = true;
        pos--;
      }
    }

    if (progress) {
      w = shortenBraid2(n, Word(std::move(w_as_vec)));
    }
  }

  std::cout << total_drop << ",";

  return w;
}

template <typename Stabilizer>
Word getReplacement(int gen) {
  if (std::is_same<Stabilizer, StabilizerDoubleSquare>::value) {
    return Word({-gen, -gen, -gen});
  }

  return Word(-gen);
}

//! Returns a pair (i, w') where i is the position where the flip was made
//! and w' is obtained from w by replacing w[i] with w[i]^{-1} or with w[i]^{-3}.
template <typename Stabilizer>
std::vector<std::pair<size_t, Word>> availableFlips(size_t n, size_t a, size_t b, Word w) {
  std::vector<std::pair<size_t, Word>> result;

  if (a > b) {
    std::swap(a, b);
  }

  // 1. Prepare initial list of words
  size_t pos = 0;
  Permutation perm1(n), perm2(n);

  for (const auto l : w) {
    const auto ind = std::abs(l);

    auto strand1 = perm2[ind - 1];
    auto strand2 = perm2[ind];

    perm1.change(strand1, strand2);
    perm2.change(ind - 1, ind);

    if (strand1 > strand2) {
      std::swap(strand1, strand2);
    }

    // available flip
    if (strand1 == a && strand2 == b) {
      auto w1 = w.subword(0, pos);
      w1 *= getReplacement<Stabilizer>(l);
      w1 *= w.subword(pos + 1, w.size());

      result.push_back(std::make_pair(pos, w1));
    }

    pos++;
  }

  // 2. Apply shortenBraid2
  return parallel::map(
      result, [n](const std::pair<size_t, Word>& p) { return std::make_pair(p.first, shortenBraid2(n, p.second)); });
}

template <typename Stabilizer, typename URNG>
Word reduce(size_t n, size_t a, size_t b, const Word& w, const std::vector<Word>& betas_conjugates, URNG& g) {
  if (a > b) {
    std::swap(a, b);
  }

  std::set<Word> checked_elements;
  std::set<Word> unchecked_elements;

  cout << "Dropped Pairs: ";
  unchecked_elements.insert(dropPairs(n, a, b, w));
  cout << endl;

  auto best_length = w.length();
  size_t iterations_since_progress = 0;
  int available_attempts = 3;

  while (!unchecked_elements.empty()) {
    // 1. Pick the best candidate
    const auto w1 = *unchecked_elements.begin();
    unchecked_elements.erase(w1);
    checked_elements.insert(w1);

    auto cur_length = w1.length();

    auto t = std::time(nullptr);
    auto tm = *std::localtime(&t);

    std::cout << ".................................................." << std::endl;
    std::cout << "Current |w| = " << w1.length() << std::put_time(&tm, ",  tm = %H-%M-%S") << std::endl;

    if (best_length > cur_length) {
      best_length = cur_length;
      iterations_since_progress = 0;
    }

    // 2. Find all flips
    const auto flips = availableFlips<Stabilizer>(n, a, b, w1);

    const auto tuple_commute = parallel::bmap(
        flips, [&](const std::pair<size_t, Word>& p) { return doesCommuteWithTuple(n, p.second, betas_conjugates); });

    if (!flips.empty()) {
      cout << "Deltas: ";

      for (size_t i = 0; i < flips.size(); ++i) {
        const auto& w2 = flips[i].second;

        const auto success = tuple_commute[i];

        const auto delta = w1.length() - w2.length();

        std::cout << delta << ",";

        if (success) {
          std::cout << endl;
          return w2;
        }

        // Add new element to unchecked_elements
        if ((checked_elements.count(w2) == 0) && (unchecked_elements.count(w2) == 0)) {
          unchecked_elements.insert(w2);
        }
      }

      std::cout << endl;
    }

    // 3. We consider ourselves stuck at this point
    if (++iterations_since_progress > 200 || unchecked_elements.empty()) {
      std::cout << "Attempt failed! " << available_attempts << " attempts left." << std::endl;

      if (available_attempts-- > 0) {
        // reset using the best candidate
        Word w2 = *checked_elements.begin();
        const auto sigma_w2 = coloredburau::permutation(n, w2);
        const auto sigma_w2_inv = sigma_w2.inverse();

        const int cloaking_center_degree = std::is_same<Stabilizer, StabilizerDoubleSquare>::value ? 4 : 2;

        for (size_t i = 0; i < 3; ++i) {
          w2 = coloredburau::generateCloakingElement(
                   n, sigma_w2[a], sigma_w2[b], 30, 50, sigma_w2_inv, cloaking_center_degree, g)
               * w2;
        }

        for (size_t i = 0; i < 3; ++i) {
          w2 = w2 * coloredburau::generateCloakingElement(n, a, b, 30, 50, sigma_w2, cloaking_center_degree, g);
        }

        // Update sets
        unchecked_elements.clear();
        checked_elements.clear();
        iterations_since_progress = 0;
        cout << "Dropped Pairs: ";
        // w2 = dropPairs(n, a, b, shortenBraid2(n, w2));
        w2 = dropPairs(n, a, b, shortBraidForm(n, w2));
        best_length = w2.length();
        unchecked_elements.insert(w2);
        cout << endl;
      } else {
        throw std::logic_error("Time expired.");
      }
    }
  }

  // Should not happen
  std::cout << "What?" << std::endl;
  throw std::logic_error("What?");
}

//! Returns true if w is written in the generators of the subgroup <b_1,...,b_{n/2 - 1}> of B_n.
bool belongsToLn(size_t n, const Word& w) {
  const auto m = occurrences(w);

  for (size_t i = n / 2; i < n; ++i) {
    if (m.count(i) != 0) {
      return false;
    }
  }

  return true;
}

//! Checks if w commutes with the subgroup <b_1,...,b_{n/2 - 1}> of B_n.
bool doesCommuteWithLn(size_t n, const Word& w) {
  std::vector<Word> l_n_genreators;
  l_n_genreators.reserve(n / 2);

  for (int i = 1; i < n / 2; ++i) {
    l_n_genreators.push_back(Word(i));
  }

  return doesCommuteWithTuple(n, w, l_n_genreators);
}

//! Returns true if w is written in the generators of the subgroup <b_{n/2 + 1},...,b_{n-1}> of B_n.
bool belongsToUn(size_t n, const Word& w) {
  const auto m = occurrences(w);

  for (size_t i = 1; i < n / 2 + 1; ++i) {
    if (m.count(i) != 0) {
      return false;
    }
  }

  return true;
}

//! Checks if w commutes with the subgroup <b_{n/2 + 1},...,b_{n-1}> of B_n.
bool doesCommuteWithUn(size_t n, const Word& w) {
  std::vector<Word> u_n_genreators;
  u_n_genreators.reserve(n / 2);

  for (int i = n / 2 + 1; i < n; ++i) {
    u_n_genreators.push_back(Word(i));
  }

  return doesCommuteWithTuple(n, w, u_n_genreators);
}

template <typename T, typename Stabilizer, typename URNG>
bool attack(const ProtocolInstance<T>& instance, URNG& g) {
  const auto& parameters = instance.parameters();

  const auto n = parameters.n();

  const auto& t_values = parameters.tValues();

  const auto& pub_a = instance.alicePublicKey();
  const auto& priv_a = instance.alicePrivateKey();

  const auto& pub_b = instance.bobPublicKey();
  const auto& priv_b = instance.bobPrivateKey();

  const auto& betas = instance.betasConjugates();

  const auto sigma_b = coloredburau::permutation(n, pub_b);

  const auto a = parameters.a();
  const auto b = parameters.b();

  // 1. Recover the conjugator z
  Word recovered_z;
  std::vector<Word> new_betas;

  std::cout << "Recovering conjugator c and betas from betas conjugates..." << std::endl;
  std::tie(recovered_z, new_betas) = recoverZ(n, betas);
  std::cout << "c and betas recovered" << std::endl;

  // DIAGNOSTICS
  std::cout << "|c^-1 * z| = " << shortenBraid2(n, -recovered_z * instance.z()).length() << std::endl;

  // recovered betas in U_n generators
  std::cout << "c^-1 * z * beta_i * z^-1 * c in U_n generators: ";
  const auto new_betas_in_U_n = parallel::bmap(new_betas, [n](const Word& b) { return belongsToUn(n, b); });

  for (const auto b : new_betas_in_U_n) {
    std::cout << (b ? "+" : "X");
  }
  std::cout << std::endl;

  // recovered betas commute with L_n
  std::cout << "c^-1 * z * beta_i * z^-1 * c commute with L_n: ";
  const auto new_betas_commute_with_L_n =
      parallel::bmap(new_betas, [n](const Word& b) { return doesCommuteWithLn(n, b); });

  for (const auto b : new_betas_commute_with_L_n) {
    std::cout << (b ? "+" : "X");
  }
  std::cout << std::endl;

  std::cout << "c^-1 * z * a * z^-1 * c in L_n generators: " << std::boolalpha
            << belongsToLn(n, shortenBraid2(n, -recovered_z * instance.alicePrivateKey() * recovered_z)) << std::endl;

  std::cout << "c^-1 * z * a * z^-1 * c commutes with U_n: " << std::boolalpha
            << doesCommuteWithUn(n, -recovered_z * instance.alicePrivateKey() * recovered_z) << std::endl;

  const auto sigma_z = coloredburau::permutation(n, recovered_z);

  //  std::cout << "New betas:" << std::endl;
  //  for (const auto& beta : new_betas) {
  //    std::cout << beta << std::endl;
  //  }

  // 2. The attack
  Word result;
  try {
    result = reduce<Stabilizer>(
        n, sigma_z[sigma_b[a]], sigma_z[sigma_b[b]], (-recovered_z) * pub_a * recovered_z, new_betas, g);
  } catch (const std::logic_error&) {
    return false;
  }

  result = recovered_z * result * -recovered_z;

  // Diagnostics: Check if the result acts as Alice's private key
  if (coloredburau::project(priv_b * result, t_values) != coloredburau::project(priv_b * pub_a, t_values)) {
    std::cout << "Restored private key doesn't act as Alice's private key" << std::endl;
    return false;
  }

  // Check if we actually found the original key
  if (areEqualBraids(n, priv_a, result)) {
    std::cout << "Recovered the original private key" << std::endl;
  }

  return true;
}

using GF32 = finitefield::FieldElement<finitefield::IdealGeneratedByPolynomial<finitefield::ZZ<2>, 1, 0, 1, 0, 0, 1>>;
using GF256 =
    finitefield::FieldElement<finitefield::IdealGeneratedByPolynomial<finitefield::ZZ<2>, 1, 1, 0, 1, 1, 0, 0, 0, 1>>;

template <typename Stabilizer = StabilizerSquare>
Protocol<GF32, GarsideDehornoyObfuscator, StochasticRewriteObfuscator, Stabilizer>
getProtocolFor128BitsSecurity(size_t seed) {
  const size_t n = 16;
  const size_t L = 15;
  const size_t r = 32;
  const size_t bob_private_key_length = 22;

  const auto random_parameters = randomParameters<GF32, Stabilizer>(n, seed)
                                     .r(r)
                                     .zMinLength(180)
                                     .zMaxLength(250)
                                     .alphaMinLength(300)
                                     .alphaMaxLength(400)
                                     .betaMinLength(50)
                                     .betaMaxLength(100)
                                     .bobPrivateKeyMinLength(bob_private_key_length)
                                     .bobPrivateKeyMaxLength(bob_private_key_length);

  const auto stabilizer = Stabilizer(random_parameters, L);

  return getProtocol(
      random_parameters, GarsideDehornoyObfuscator(), getStochasticRewriteObfuscator(n, seed), stabilizer);
}

template <typename Stabilizer = StabilizerSquare>
Protocol<GF256, GarsideDehornoyObfuscator, StochasticRewriteObfuscator, Stabilizer>
getProtocolFor256BitsSecurity(size_t seed) {
  const size_t n = 16;
  const size_t L = 30;
  const size_t r = 32;
  const size_t bob_private_key_length = 43;

  const auto random_parameters = randomParameters<GF256, Stabilizer>(n, seed)
                                     .r(r)
                                     .zMinLength(300)
                                     .zMaxLength(400)
                                     .alphaMinLength(300)
                                     .alphaMaxLength(400)
                                     .betaMinLength(100)
                                     .betaMaxLength(200)
                                     .bobPrivateKeyMinLength(bob_private_key_length)
                                     .bobPrivateKeyMaxLength(bob_private_key_length);

  const auto stabilizer = Stabilizer(random_parameters, L);

  return getProtocol(
      random_parameters, GarsideDehornoyObfuscator(), getStochasticRewriteObfuscator(n, seed), stabilizer);
}

Protocol<GF32, GarsideDehornoyObfuscator, StochasticRewriteObfuscator, StabilizerManyShort>
getProtocolFor128BitsSecurityMultipleCloaking(size_t seed) {
  const size_t n = 16;
  const size_t r = 32;
  const size_t bob_private_key_length = 22;

  const size_t cloak_min_length = 30;
  const size_t cloak_max_length = 50;
  const size_t cloaking_number = 30;

  const auto random_parameters = randomParameters<GF32, StabilizerManyShort>(n, seed)
                                     .r(r)
                                     .zMinLength(180)
                                     .zMaxLength(250)
                                     .alphaMinLength(300)
                                     .alphaMaxLength(400)
                                     .betaMinLength(50)
                                     .betaMaxLength(100)
                                     .bobPrivateKeyMinLength(bob_private_key_length)
                                     .bobPrivateKeyMaxLength(bob_private_key_length);

  const auto stabilizer = StabilizerManyShort(random_parameters, cloak_min_length, cloak_max_length, cloaking_number);

  return getProtocol(
      random_parameters, GarsideDehornoyObfuscator(), getStochasticRewriteObfuscator(n, seed), stabilizer);
}

Protocol<GF256, GarsideDehornoyObfuscator, StochasticRewriteObfuscator, StabilizerManyShort>
getProtocolFor256BitsSecurityMultipleCloaking(size_t seed) {
  const size_t n = 16;
  const size_t r = 32;
  const size_t bob_private_key_length = 43;

  const size_t cloak_min_length = 30;
  const size_t cloak_max_length = 50;
  const size_t cloaking_number = 60;

  const auto random_parameters = randomParameters<GF256, StabilizerManyShort>(n, seed)
                                     .r(r)
                                     .zMinLength(300)
                                     .zMaxLength(400)
                                     .alphaMinLength(300)
                                     .alphaMaxLength(400)
                                     .betaMinLength(100)
                                     .betaMaxLength(200)
                                     .bobPrivateKeyMinLength(bob_private_key_length)
                                     .bobPrivateKeyMaxLength(bob_private_key_length);

  const auto stabilizer = StabilizerManyShort(random_parameters, cloak_min_length, cloak_max_length, cloaking_number);

  return getProtocol(
      random_parameters, GarsideDehornoyObfuscator(), getStochasticRewriteObfuscator(n, seed), stabilizer);
}
} // namespace kayawood
} // namespace crag

#endif // CRAG_KAYAWOOD_ATTACK_H
