/*
 * permutation16.h
 *
 *  Created on: Apr 28, 2013
 *      Author: dpantele
 */

#ifndef PERMUTATION16_H_
#define PERMUTATION16_H_

#include <cstdint>
#include <initializer_list>
#include <random>

#include "common.h"

namespace crag {
//! This class represents the permutation of 16 elements
class Permutation16 {
  public:
    CONSTEXPR_OR_CONST static size_t RANK = 16;
    //!Trivial permutation
    CONSTEXPR Permutation16()
      : permutation_(TRIVIAL_PERMUTATION_REPRESENTATION)
    { }

    //! Define permutation with the images of 0, 1, 2... (note that elements starts from 0)
    template <class ForwardIterator>
    Permutation16(ForwardIterator begin, ForwardIterator end)
      : permutation_(0)
    {
      size_t shift = 0;
      for (auto image = begin; image != end; ++image) {
        assert(0 <= *image && *image < 16);
        assert(shift < 64);
        permutation_ |= (static_cast<uint64_t>(*image) << shift);
        shift += 4;
      }
      if (shift != 64) {
        permutation_ |= (TRIVIAL_PERMUTATION_REPRESENTATION & ((0xffffffffffffffffull) << shift));
      }
    }

    Permutation16(std::initializer_list<uint64_t> permutation)
      : Permutation16(permutation.begin(), permutation.end())
    { }

    //!Compare two permutations
    bool operator==(const Permutation16& other) const {
      return permutation_ == other.permutation_;
    }

    //Why always follow the pattern?
    bool operator!=(const Permutation16& other) const {
      return permutation_ != other.permutation_;
    }

    CONSTEXPR static size_t size() {
      return RANK;
    }

    //! Get the image of i
    uint64_t get_image(size_t i) const {
      i = i << 2; //multiple by 4
      return (permutation_ >> i) & 0xfull;
    }

    //! alias for get_image
    uint64_t operator[](size_t i) const {
      return get_image(i);
    }

    //! Set the image of element
    Permutation16& set_image(size_t element, uint64_t image) {
      assert(element < 16);
      assert(image < 16);
      element = element << 2;
      permutation_ &= ~(0xfull << element);
      permutation_ |= (image << element);

      return *this;
    }

    //! Compose this with other, this[i] = other[this[i]]
    Permutation16& compose_with(const Permutation16& other) {
      uint64_t current_permutation = permutation_;
      permutation_ = 0;

      for (size_t shift = 0; shift < 64; shift += 4) {
        permutation_ |= (other.get_image(0xfull & current_permutation) << shift);
        current_permutation = current_permutation >> 4;
      }

      return *this;
    }

    Permutation16& operator*=(const Permutation16& other) {
      return compose_with(other);
    }

    Permutation16 operator*(const Permutation16& other) const {
      return Permutation16(*this).compose_with(other);
    }


    //! Inverse of this
    Permutation16 inverse() const {
      uint64_t permutation_inverse_ = 0;

      for (uint64_t element = 0; element < 16; ++element) {
        permutation_inverse_ |= (element << (this->get_image(element) << 2));
      }

      return Permutation16(permutation_inverse_);
    }

    //!Swap the images of two elements
    Permutation16& swap_images(size_t i, size_t j) {
      if (i == j) {
        return *this;
      }

      i = i << 2;
      j = j << 2;
      permutation_ ^= ((permutation_ >> j) & 0xfull) << i;
      permutation_ ^= ((permutation_ >> i) & 0xfull) << j;
      permutation_ ^= ((permutation_ >> j) & 0xfull) << i;

      return *this;
    }

    template <class RandomEngine>
    static Permutation16 random(size_t max_element, RandomEngine& engine) {
      assert(max_element <= 16);
      Permutation16 result;

      for(size_t i = 0; i < max_element; ++i) {
        size_t swap_position = i + engine() % (max_element + 1 - i);
        result.swap_images(i, swap_position);
      }

      return result;
    }

    static Permutation16 random(size_t max_element) {
      static std::mt19937 default_engine(rand() != 0 ? rand() : 17);
      return random(max_element, default_engine);
    }

    size_t permutation_hash() const {
      CONSTEXPR_OR_CONST static std::hash<uint64_t> uint_hasher_ = std::hash<uint64_t>();
      return uint_hasher_(permutation_);
    }

  private:
    const static uint64_t TRIVIAL_PERMUTATION_REPRESENTATION = 0xfedcba9876543210ull;
    uint64_t permutation_;

    Permutation16(uint64_t permutation)
      : permutation_(permutation)
    { }
};

}

namespace std {

template<>
struct hash<crag::Permutation16> {
  public:
    size_t operator()(const crag::Permutation16& permutation) const {
      return permutation.permutation_hash();
    }
};

}
#endif /* PERMUTATION16_H_ */
