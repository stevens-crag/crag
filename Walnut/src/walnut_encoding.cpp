#include "walnut_encoding.h"

#include <random>

namespace crag {
namespace walnut {

Word getFreePureBraidSubgroupGen(size_t n, size_t i) {
  if (i + 1 > n) {
    throw std::invalid_argument("Index i is too large.");
  }

  Word prefix;

  for (size_t j = n - 1; j > i; --j) {
    prefix.push_back(j);
  }

  return prefix * (Word(i) ^ 2) * -prefix;
}

bool areValidEncodingIndices(size_t n, const std::vector<size_t>& encoding_indices) {
  if (encoding_indices.size() != 4) {
    return false;
  }

  std::set<size_t> indices(encoding_indices.begin(), encoding_indices.end());

  if (indices.size() != encoding_indices.size()) {
    return false;
  }

  if ((*indices.begin() == 0) || (*indices.rbegin() >= n)) {
    return false;
  }

  return true;
}

//! Transforms lower 4 bits into generator and exponent for Walnut.
static std::pair<int, int> toGeneratorAndExponent(uint8_t b) {
  return std::make_pair((b & 0b11) + 1, ((b >> 2) & 0b11) + 1);
}

//! Encodes lower 4 bits
static Word encodeLower(uint8_t b) {
  b = b & 0b1111;

  int gen, exp;
  std::tie(gen, exp) = toGeneratorAndExponent(b);

  return Word(gen) ^ exp;
}

Word encode(uint8_t b) {
  return encodeLower(b >> 4) * encodeLower(b);
}

Word encode(const msg_hash_t& msg_hash) {
  Word result;

  for (const auto b : msg_hash) {
    result *= encode(b);
  }

  return result;
}

Word encode(size_t n, const msg_hash_t& msg_hash, const std::vector<size_t>& encoding_indices) {
  if (n < 5) {
    throw std::invalid_argument("Require n >= 5.");
  }

  const std::vector<Word> random_gens = {
      getFreePureBraidSubgroupGen(n, encoding_indices[0]),
      getFreePureBraidSubgroupGen(n, encoding_indices[1]),
      getFreePureBraidSubgroupGen(n, encoding_indices[2]),
      getFreePureBraidSubgroupGen(n, encoding_indices[3]),
  };

  const auto w = encode(msg_hash);

  return w.replaceGenerators(random_gens);
}

msg_hash_t randomMessageHash(size_t n, size_t seed) {
  assert(n % 8 == 0);

  std::mt19937_64 g(seed);
  return randomMessageHash(n, g);
}

DefaultEncoder::DefaultEncoder(size_t n, std::vector<size_t> encoding_indices, size_t hash_size)
    : encoding_indices_(std::move(encoding_indices))
    , free_sugroup_gens_(getFreeSubgroupGens_(n, encoding_indices_))
    , hash_size_(hash_size) {
  if (hash_size % 8 != 0) {
    throw std::invalid_argument("Invalid hash size.");
  }
}

const std::vector<size_t>& DefaultEncoder::encodingIndices() const {
  return encoding_indices_;
}

Word DefaultEncoder::operator()(const msg_hash_t& msg_hash) const {
  if (hash_size_ / 8 != msg_hash.size()) {
    throw std::invalid_argument("Invalid message hash size.");
  }

  const auto w = encode(msg_hash);
  return w.replaceGenerators(free_sugroup_gens_);
}

size_t DefaultEncoder::hashSize() const {
  return hash_size_;
}

std::vector<Word> DefaultEncoder::getFreeSubgroupGens_(size_t n, const std::vector<size_t>& encoding_indices) {
  if (!areValidEncodingIndices(n, encoding_indices)) {
    throw std::invalid_argument("Invalid encoding indices.");
  }

  return {
      getFreePureBraidSubgroupGen(n, encoding_indices[0]),
      getFreePureBraidSubgroupGen(n, encoding_indices[1]),
      getFreePureBraidSubgroupGen(n, encoding_indices[2]),
      getFreePureBraidSubgroupGen(n, encoding_indices[3]),
  };
}

DefaultEncoder randomEncoder(size_t n, size_t hash_size, size_t seed) {
  std::mt19937_64 g(seed);
  return randomEncoder(n, hash_size, g);
}
} // namespace walnut
} // namespace crag
