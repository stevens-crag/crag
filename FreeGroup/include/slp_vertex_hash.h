/**
 * \file slp_vertex_hash.h
 * \brief Automorphisms from SLP result to some small group used as hashed and algorithms with them.
 */

#ifndef SLP_VERTEX_HASH_H_
#define SLP_VERTEX_HASH_H_

#include <cassert>
#include <new>
#include <memory>

#include "gmpxx.h"

typedef mpz_class LongInteger;

#include "Permutation.h"
#include "slp_vertex.h"
#include "slp_reduce.h"

namespace crag {
namespace slp {

//! Basic class to calculate hashes. Hash types are listed as template parameters.
//template <class... Hashers> class TVertexHash {
//  public:
//    /**
//     * Calculate hash of some terminal vertex using its singed id, usually
//     * just conversion of value to Vertex::SignedVertexId. All hashers should implement such constructor.
//     */
//    TVertexHash(Vertex::VertexSignedId terminal_id);
//
//    //! Hash of empty vertex. Also should be implemented by every hasher.
//    TVertexHash();
//
//    //! Copy constructor. Should be implemented by each hasher.
//    TVertexHash(const TVertexHash& other);
//
//    //! Calculate the hash of concatenated words. Each hasher should implement is as concatenate_with method.
//    TVertexHash& operator*=(const TVertexHash& other);
//
//    //! Hash comparison. Each hasher should implement it as method is_equal_to
//    bool operator==(const TVertexHash& other);
//
//    //! Return the hash of inverted word.
//    TVertexHash inverse() const;
//
//    //! Replace with the hash of inverted word and return *this. Each hasher should implements it.
//    TVertexHash& inverse_inplace();
//};
//
//Here are details of realisation, no documentation
//template <class TFirstHasher, class... TOtherHashers>
//class TVertexHash<TFirstHasher, TOtherHashers...> : public TFirstHasher, public TVertexHash<TOtherHashers...> {
//    typedef TFirstHasher FirstHasher;
//    typedef TVertexHash<TOtherHashers...> OtherHasher;
//  public:
//
//    TVertexHash(Vertex::VertexSignedId terminal_id)
//      : FirstHasher(terminal_id)
//      , OtherHasher(terminal_id)
//    { }
//
//    TVertexHash()
//    { }
//
//    TVertexHash(const TVertexHash& other)
//      : FirstHasher(other)
//      , OtherHasher(other)
//    { }
//
//    TVertexHash& operator*=(const TVertexHash& other) {
//      FirstHasher::concatenate_with(other);
//      OtherHasher::operator*=(other);
//
//      return *this;
//    }
//
//    //I do not add move support right now, since for current hashers it makes no sense
//
//    bool operator==(const TVertexHash& other) {
//      return FirstHasher::is_equal_to(other) && OtherHasher::operator==(other);
//    }
//
//    TVertexHash inverse() const {
//      TVertexHash copy(*this);
//      return copy.inverse_inplace();
//    }
//
//    TVertexHash& inverse_inplace() {
//      FirstHasher::inverse_inplace();
//      OtherHasher::inverse_inplace();
//
//      return *this;
//    }
//
//    typedef std::unordered_map<Vertex, TVertexHash> Cache;
//
//
//};
//
//template <>
//class TVertexHash<> {
//  public:
//    TVertexHash& operator*=(const TVertexHash& other) {
//      return *this;
//    }
//
//    TVertexHash& inverse_inplace() {
//      return *this;
//    }
//
//    TVertexHash inverse() const {
//      return *this;
//    }
//
//    TVertexHash(Vertex::VertexSignedId terminal) {}
//    TVertexHash() {}
//
//    bool operator==(const TVertexHash& other) {
//      return true;
//    }
//};

union MaxAlign
{
    int                 i     ;
    long                l     ;
    long long           ll    ;
    long double         ld    ;
    double              d     ;
    void*               p     ;
    void (*             pf)() ;
    MaxAlign*           ps    ;
} ;

template <class... Hashers> class TVertexHash;

template <class... Hashers> class TVertexHash;
template <class TFirstHasher, class... TOtherHashers>
class TVertexHash<TFirstHasher, TOtherHashers...> : public TVertexHash<TOtherHashers...> {
    typedef TFirstHasher FirstHasher;
    typedef TVertexHash<TOtherHashers...> OtherHasher;
  public:
    TVertexHash()
      : OtherHasher()
      , first_hasher_(nullptr)
    {}

    TVertexHash(Vertex::VertexSignedId terminal_symbol)
      : OtherHasher(terminal_symbol)
      , first_hasher_(nullptr)
    { }

    TVertexHash(std::shared_ptr<TVertexHash> left, std::shared_ptr<TVertexHash> right)
      : OtherHasher(left, right)
      , first_hasher_(nullptr)
    { }

    TVertexHash(const TVertexHash& other)
      : OtherHasher(other)
      , first_hasher_(nullptr)
    {
      if (other.first_hasher_) {
        first_hasher_ = reinterpreter_cast<FirstHasher*>(new (first_hasher_storage_.o) FirstHasher(*(other.first_hasher_)));
      }
    }

    ~TVertexHashImpl() {
      first_hasher_->~FirstHasher();
    }

    const FirstHasher& calculate_first(Vertex::VertexSignedId terminal_id) const {
      if (!first_hasher_) {
        first_hasher_ = reinterpreter_cast<FirstHasher*>(new (first_hasher_storage_.o) FirstHasher(terminal_id));
      }
      return *first_hasher_;
    }

    const FirstHasher& calculate_first(const FirstHasher& left, const FirstHasher& right) const {
      if (!first_hasher_) {
        first_hasher_ = reinterpreter_cast<FirstHasher*>(new (first_hasher_storage_.o) FirstHasher(left));
        *first_hasher_ *= right;
      }
      return *first_hasher_;
    }

    const FirstHasher& calculate_first(const TVertexHash& this_container) const {
      if (!first_hasher_) {
        first_hasher_ = reinterpreter_cast<FirstHasher*>(new (first_hasher_storage_.o) FirstHasher(terminal_id));
      }
      return *first_hasher_;
    }

    bool is_equal_to(const TVertexHash& this_container, const TVertexHash& other_container) const;

    TVertexHash inverse() const {
      TVertexHash copy(*this);
      return copy.inverse_inplace();
    }

    TVertexHash& inverse_inplace() {
      FirstHasher::inverse_inplace();
      OtherHasher::inverse_inplace();

      return *this;
    }

    typedef std::unordered_map<Vertex, TVertexHash> Cache;

  private:
    union
    {
        MaxAlign dummyForAlignment;
        unsigned char o[sizeof(FirstHasher)];
    } first_hasher_storage_;
    mutable FirstHasher* first_hasher_;
};

template <>
class TVertexHash<> {
  public:
    TVertexHash& operator*=(const TVertexHash& other) {
      return *this;
    }

    TVertexHash& inverse_inplace() {
      return *this;
    }

    TVertexHash inverse() const {
      return *this;
    }

    TVertexHash(Vertex::VertexSignedId terminal) {}
    TVertexHash() {}

    bool operator==(const TVertexHash& other) {
      return true;
    }
  protected:
    std::shared_ptr<TVertexHash> left_hash_or_negate_;
    std::shared_ptr<TVertexHash> right_hash_or_null_;
    Vertex::VertexSignedId terminal_symbol_;
};

template <class... Hashers> class TVertexHash {
    template <class...> friend TVertexHashImpl;
  private:
    TVertexHashImpl<Hashers> impl;
    Vertex::VertexSignedId terminal_symbol_;
};



namespace hashers {

//! Calculate the power of each terminal assuming considering the group as commutative.
template <size_t RANK>
class PowerCountHash {
  private:
    int64_t terminal_power_[RANK] = {0};
  public:
    PowerCountHash() {}
    PowerCountHash(Vertex::VertexSignedId terminal_id) {
      assert(terminal_id < RANK && -terminal_id < RANK && terminal_id != 0);
      if (terminal_id > 0) {
        terminal_power_[terminal_id - 1] = 1;
      } else if (terminal_id < 0){
        terminal_power_[-terminal_id - 1] = -1;
      }
    }

    PowerCountHash(const PowerCountHash& other) {
      for (size_t i = 0; i < RANK; ++i) {
        terminal_power_[i] = other.terminal_power_[i];
      }
    }

    void concatenate_with(const PowerCountHash& other) {
      for (size_t i = 0; i < RANK; ++i) {
        terminal_power_[i] += other.terminal_power_[i];
      }
    }

    void inverse_inplace() {
      for (size_t i = 0; i < RANK; ++i) {
        terminal_power_[i] *= -1;
      }
    }

    bool is_equal_to(const PowerCountHash& other) {
      for (size_t i = 0; i < RANK; ++i) {
        if (terminal_power_[i] != other.terminal_power_[i]) {
          return false;
        }
      }
      return true;
    }
};

class SinglePowerHash {
  private:
    int64_t terminals_power_;
  public:
    SinglePowerHash()
      : terminals_power_(0)
    {}
    SinglePowerHash(Vertex::VertexSignedId terminal_id)
      : terminals_power_(terminal_id > 0 ? 1 : (terminal_id < 0 ? -1 : 0))
    { }

    SinglePowerHash(const SinglePowerHash& other) {
      terminals_power_ = other.terminals_power_;
    }

    void concatenate_with(const SinglePowerHash& other) {
      terminals_power_ += other.terminals_power_;
    }

    void inverse_inplace() {
      terminals_power_ *= -1;
    }

    bool is_equal_to(const SinglePowerHash& other) {
      return terminals_power_ == other.terminals_power_;
    }
};


//! Assign some random permutations to each of the terminals and work in the group of substitutions.
class PermutationHash {
  private:
    constexpr static size_t PERMUTATION_RANK = 16;
    Permutation permutation_;
    static Permutation GetTerminalPermutation(Vertex::VertexSignedId terminal_id) {
      static std::vector<Permutation> permutations = {
        Permutation(PERMUTATION_RANK), //for null terminal
        Permutation(std::vector<int>({11, 4, 5, 13, 15, 0, 12, 8, 3, 1, 6, 14, 9, 7, 2, 10})), //permutation of the maximal order in S16
        Permutation(std::vector<int>({6, 14, 0, 4, 13, 7, 11, 12, 1, 10, 15, 9, 5, 8, 2, 3})), //combined with the previous, it can give the whole group
      };

      size_t terminal = terminal_id < 0 ? -terminal_id : terminal_id;

      while (permutations.size() <= terminal) {
        permutations.push_back(Permutation::random(PERMUTATION_RANK));
      }

      if (terminal_id >= 0) {
        return permutations[terminal_id];
      } else {
        return permutations[-terminal_id].inverse();
      }
    }

  public:
    PermutationHash(Vertex::VertexSignedId terminal_id)
      : permutation_(GetTerminalPermutation(terminal_id))
    { }

    PermutationHash()
      : permutation_(GetTerminalPermutation(0))
    { }

    PermutationHash(const PermutationHash& other)
      : permutation_(other.permutation_)
    { }

    void concatenate_with(const PermutationHash& other) {
      permutation_ *= other.permutation_;
    }

    void inverse_inplace() {
      permutation_ = permutation_.inverse();
    }

    bool is_equal_to(const PermutationHash& other) {
      return permutation_ == other.permutation_;
    }
};

} //namespace hashers

//! Here some algorithms which use hashed are implemented. Made for the convenience when using templates.
template <class... Hashers>
class TVertexHashAlgorithms {
  public:
    typedef TVertexHash<Hashers...> VertexHash; //!< The type of basic hash.

    /**
     * All algorithms here can take a pointer to the object of this type to store calculate hashes.
     */
    typedef std::unordered_map<Vertex, VertexHash> Cache;

    //! Calculate the hash of the word produced by root, starting from begin to end.
    static VertexHash get_subvertex_hash(const Vertex& root, const LongInteger& begin, const LongInteger& end, Cache* cache) {
      assert(cache);
      assert(begin >= 0 && end <= root.length());

      static VertexHash Null;
      if (begin >= root.length() || end < 0 || end <= begin) {
        return Null;
      }
      if (root.height() == 1) {
        return cache->insert(
            std::make_pair(
                root,
                VertexHash(root.vertex_id())
            )
        ).first->second;
      } else {
        if (begin == 0 && end == root.length()) {
          auto cache_item = cache->find(root);

          if (cache_item != cache->end()) {
            return cache_item->second;
          } else {
            cache_item = cache->find(root.negate());
            if (cache_item != cache->end()) {
              return cache_item->second.inverse();
            }

            const VertexHash& result = cache->insert(
                std::make_pair(
                    root,
                    VertexHash(get_subvertex_hash(root.left_child(), begin, root.split_point(), cache)) *=
                               get_subvertex_hash(root.right_child(), 0, end - root.split_point(), cache)
                )).first->second;

            return result;
          }
        }
        if (root.split_point() >= end) {
          return get_subvertex_hash(root.left_child(), begin, end, cache);
        } else if (root.split_point() <= begin) {
          return get_subvertex_hash(root.right_child(), begin - root.split_point(), end - root.split_point(), cache);
        } else {
          return VertexHash(get_subvertex_hash(root.left_child(), begin, root.split_point(), cache)) *=
                            get_subvertex_hash(root.right_child(), 0, end - root.split_point(), cache);
        }
      }
    }

    //! Call get_subvertex_hash(root, begin, end, temporary cache)
    static VertexHash get_subvertex_hash(const Vertex& root, const LongInteger& begin, const LongInteger& end) {
      Cache cache;
      return get_subvertex_hash(root, begin, end, &cache);
    }

    //! Get the longest common prefix using binary search and hashes
    static LongInteger get_longest_common_prefix(
        const Vertex& first,
        const Vertex& second,
        Cache* calculated_hashes
    ) {

      if (first.length() > second.length()) {
        return get_longest_common_prefix(second, first, calculated_hashes);
      }

      LongInteger begin = 0;
      LongInteger end = first.length();

      while (begin < end) {
        LongInteger split = (begin + end) / 2 + 1;
        if (get_subvertex_hash(first, begin, split, calculated_hashes) == get_subvertex_hash(second, begin, split, calculated_hashes)) {
          begin = split;
        } else {
          end = split - 1;
        }
      }

      return begin;
    }

    //! Call get_longest_common_prefix(first, second, temporary cache)
    static LongInteger get_longest_common_prefix(
            const Vertex& first,
            const Vertex& second) {
      Cache cache;
      return get_longest_common_prefix(first, second, &cache);
    }


    static LongInteger get_cancellation_length(
        const Vertex& vertex,
        Cache* calculated_hashes) {
      return get_longest_common_prefix(
          vertex.left_child().negate(),
          vertex.right_child(),
          calculated_hashes
      );
    }

    //! Call get_cancellation_length(vertex, temporary cache)
    static LongInteger get_cancellation_length(const Vertex& vertex) {
      Cache cache;
      return get_cancellation_length(vertex, &cache);
    }

    static Vertex reduce(
        const Vertex& vertex,
        Cache* calculated_hashes,
        std::unordered_map<Vertex, Vertex>* reduced_vertices
    ) {
      return base_reduce(
          vertex,
          [calculated_hashes](const Vertex& vertex) { return get_cancellation_length(vertex, calculated_hashes); },
          reduced_vertices
      );
    }

    static Vertex reduce(const Vertex& vertex) {
      Cache cache;
      std::unordered_map<Vertex, Vertex> reduced_vertices;
      return reduce(vertex, &cache, &reduced_vertices);
    }
};

} //namespace slp
} //namespace crag


#endif /* SLP_VERTEX_HASH_H_ */
