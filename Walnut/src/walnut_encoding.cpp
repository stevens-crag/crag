#include "walnut_encoding.h"

#include <random>

namespace crag {
namespace walnut {

Word getFreePureBraidSubgroupGen(size_t n, size_t i) {
  if ((i == 0) || (i + 1 > n)) {
    throw std::invalid_argument("Index i is out of range.");
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

AdvancedEncoder::AdvancedEncoder(size_t n, size_t hash_size, const std::vector<std::vector<size_t>>& generator_indices)
    : gen_matrix_(getFreeGenerators_(n, generator_indices))
    , hash_size_(hash_size) {}

Word AdvancedEncoder::operator()(const msg_hash_t& msg_hash) const {
  if (hash_size_ / 8 != msg_hash.size()) {
    throw std::invalid_argument("Invalid message hash size.");
  }

  Word result;

  const auto k = gen_matrix_.size();
  size_t i = 0;

  for (const auto b : msg_hash) {
    for (const auto idx : this->encode(b)) {
      result *= gen_matrix_[i % k][idx];
      ++i;
    }
  }

  return result;
}

size_t AdvancedEncoder::hashSize() const {
  return hash_size_;
}

std::vector<std::vector<Word>>
AdvancedEncoder::getFreeGenerators_(size_t n, const std::vector<std::vector<size_t>>& generator_indices) const {
  std::vector<std::vector<Word>> result;
  result.reserve(generator_indices.size());

  for (const auto& v : generator_indices) {
    if (v.size() != 4) {
      throw std::invalid_argument("Invalid number of indices.");
    }

    result.push_back({
        getFreePureBraidSubgroupGen(n, v[0]),
        getFreePureBraidSubgroupGen(n, v[1]),
        getFreePureBraidSubgroupGen(n, v[2]),
        getFreePureBraidSubgroupGen(n, v[3]),
    });
  }

  return result;
}

std::vector<size_t> AdvancedEncoder::encode(uint8_t b) const {
  return {
      static_cast<size_t>((b >> 6) & 0b11),
      static_cast<size_t>((b >> 4) & 0b11),
      static_cast<size_t>((b >> 2) & 0b11),
      static_cast<size_t>(b & 0b11),
  };
}

AdvancedEncoder randomAdvancedEncoder(size_t n, size_t k, size_t hash_size, size_t seed) {
  std::mt19937_64 g(seed);
  return randomAdvancedEncoder(n, k, hash_size, g);
}
} // namespace walnut
} // namespace crag
