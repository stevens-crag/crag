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
#include "permutation16.h"
#include "slp_vertex.h"
#include "slp_reduce.h"

namespace crag {
namespace slp {

//! Basic class to calculate hashes. Hash types are listed as template parameters.
template <class... Hashers> class TVertexHash {
  public:
    /**
     * Calculate hash of some terminal vertex using its singed id, usually
     * just conversion of value to Vertex::SignedVertexId. All hashers should implement such constructor.
     */
    TVertexHash(Vertex::VertexSignedId terminal_id);

    //! Hash of empty vertex. Also should be implemented by every hasher.
    TVertexHash();

    //! Copy constructor. Should be implemented by each hasher.
    TVertexHash(const TVertexHash& other);

    //! Calculate the hash of concatenated words. Each hasher should implement is as concatenate_with method.
    TVertexHash& operator*=(const TVertexHash& other);

    //! Hash comparison. Each hasher should implement it as method is_equal_to
    bool operator==(const TVertexHash& other) const;
    bool operator!=(const TVertexHash& other) const;

    //! Return the hash of inverted word.
    TVertexHash inverse() const;

    //! Replace with the hash of inverted word and return *this. Each hasher should implements it.
    TVertexHash& inverse_inplace();

    //! Return size_t value to be used as hash in std::unordered_map/set
    size_t get_std_hash() const;
};

//Here are details of implementation, no documentation
template <class TFirstHasher, class... TOtherHashers>
class TVertexHash<TFirstHasher, TOtherHashers...> : public TFirstHasher, public TVertexHash<TOtherHashers...> {
    typedef TFirstHasher FirstHasher;
    typedef TVertexHash<TOtherHashers...> OtherHasher;
  public:

    TVertexHash(Vertex::VertexSignedId terminal_id)
      : FirstHasher(terminal_id)
      , OtherHasher(terminal_id)
    { }

    TVertexHash()
    { }

    TVertexHash(const TVertexHash& other)
      : FirstHasher(other)
      , OtherHasher(other)
    { }

    TVertexHash& operator*=(const TVertexHash& other) {
      FirstHasher::concatenate_with(other);
      OtherHasher::operator*=(other);

      return *this;
    }

    //I do not add move support right now, since for current hashers it makes no sense

    bool operator==(const TVertexHash& other) const {
      return FirstHasher::is_equal_to(other) && OtherHasher::operator==(other);
    }

    bool operator!=(const TVertexHash& other) const {
      return !(*this == other);
    }

    TVertexHash inverse() const {
      TVertexHash copy(*this);
      return copy.inverse_inplace();
    }

    TVertexHash& inverse_inplace() {
      FirstHasher::inverse_inplace();
      OtherHasher::inverse_inplace();

      return *this;
    }

    size_t get_std_hash() const {
      size_t first_hash_value = FirstHasher::get_std_hash();
      return OtherHasher::get_std_hash() + 0x9e3779b9 + (first_hash_value << 6) + (first_hash_value >> 2);
    }

    typedef std::unordered_map<Vertex, TVertexHash> Cache;


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

    bool operator==(const TVertexHash& other) const {
      return true;
    }

    bool operator!=(const TVertexHash& other) const {
      return false;
    }

    size_t get_std_hash() const {
      return 0;
    }

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

    bool is_equal_to(const PowerCountHash& other) const {
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
    {}

    SinglePowerHash(const SinglePowerHash& other) {
      terminals_power_ = other.terminals_power_;
    }

    void concatenate_with(const SinglePowerHash& other) {
      terminals_power_ += other.terminals_power_;
    }

    void inverse_inplace() {
      terminals_power_ *= -1;
    }

    bool is_equal_to(const SinglePowerHash& other) const {
      return terminals_power_ == other.terminals_power_;
    }

    size_t get_std_hash() const {
      CONSTEXPR_OR_CONST static std::hash<int64_t> int64_t_hasher_ = std::hash<int64_t>();
      return int64_t_hasher_(terminals_power_);
    }
};


//! Assign some random permutations to each of the terminals and work in the group of substitutions.
template <class TPermutation, class BasePermutations>
class PermutationHashBase {
  private:
    TPermutation permutation_;
  
    static TPermutation GetTerminalPermutation(Vertex::VertexSignedId terminal_id) {
      static std::vector<TPermutation> permutations(BasePermutations::permutations());

      size_t terminal = terminal_id < 0 ? -terminal_id : terminal_id;

      while (permutations.size() <= terminal) {
        permutations.push_back(TPermutation::random(TPermutation::RANK));
      }

      if (terminal_id >= 0) {
        return permutations[terminal_id];
      } else {
        return permutations[-terminal_id].inverse();
      }
    }

  public:
    PermutationHashBase(Vertex::VertexSignedId terminal_id)
      : permutation_(GetTerminalPermutation(terminal_id))
    { }

    PermutationHashBase()
      : permutation_(GetTerminalPermutation(0))
    { }

    PermutationHashBase(const PermutationHashBase& other)
      : permutation_(other.permutation_)
    { }

    void concatenate_with(const PermutationHashBase& other) {
      permutation_ *= other.permutation_;
    }

    void inverse_inplace() {
      permutation_ = permutation_.inverse();
    }

    bool is_equal_to(const PermutationHashBase& other) const {
      return permutation_ == other.permutation_;
    }

    size_t get_std_hash() const {
      CONSTEXPR_OR_CONST static std::hash<TPermutation> permutation_hasher_ = std::hash<TPermutation>();
      return permutation_hasher_(permutation_);
    }

};

template <class TPermutation>
class DefaultBasePermutations {
  public:
    static const std::vector<TPermutation>& permutations() {
      static std::vector<TPermutation> permutations = {
        TPermutation(), //for null terminal
        TPermutation({11, 4, 5, 13, 15, 0, 12, 8, 3, 1, 6, 14, 9, 7, 2, 10}), //permutation of the maximal order in S16
        TPermutation({6, 14, 0, 4, 13, 7, 11, 12, 1, 10, 15, 9, 5, 8, 2, 3}), //combined with the previous, it can give the whole group
      };
      return permutations;
    }
};


//! The hash with the given default permutations.
template <class TPermutation>
using PermutationHash = PermutationHashBase<TPermutation, DefaultBasePermutations<TPermutation>>;

class ImageLengthHash {
  private:
    LongInteger length_;
  public:
    ImageLengthHash()
      : length_(0)
    {}
    ImageLengthHash(Vertex::VertexSignedId terminal_id)
      : length_(terminal_id != 0 ? 1 : 0)
    {}

    ImageLengthHash(const ImageLengthHash& other)
      : length_(other.length_)
    {}

    void concatenate_with(const ImageLengthHash& other) {
      length_ += other.length_;
    }

    void inverse_inplace()
    {}

    bool is_equal_to(const ImageLengthHash& other) const {
      return length_ == other.length_;
    }

    size_t get_std_hash() const {
      CONSTEXPR_OR_CONST static std::hash<LongInteger> long_integer_hasher_ = std::hash<LongInteger>();
      return long_integer_hasher_(length_);
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

    /**
     * Type of cache of vertices representing given hashes. It is used to remove duplicates.
     */
    typedef std::unordered_map<VertexHash, Vertex> HashRepresentativesCache;

    //! Calculate the hash of the word produced by root, starting from begin to end, filling hash cache for all the subvertices.
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

    //! Calculate the hash of the word produced by root, filling hash cache for all the subvertices
    static VertexHash get_vertex_hash(const Vertex& root, Cache* cache) {
      return get_subvertex_hash(root, 0, root.length(), cache);
    }

    //! Calculate the hash of the word produced by root using temporary cache.
    static VertexHash get_vertex_hash(const Vertex& root) {
      return get_subvertex_hash(root, 0, root.length());
    }

    //! Get the longest common prefix using binary search and hashes
    static LongInteger get_longest_common_prefix(
        const Vertex& first,
        const Vertex& second,
        LongInteger begin, //already verified left boundary
        LongInteger end, //already verified right boundary
        Cache* calculated_hashes
    ) {

      if (first.length() > second.length()) {
        return get_longest_common_prefix(second, first, std::move(begin), std::move(end), calculated_hashes);
      }

      if (end == -1 || end > first.length()) {
        end = first.length();
      }

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

    static LongInteger get_longest_common_prefix(
        const Vertex& first,
        const Vertex& second,
        Cache* calculated_hashes
    ) {
      return get_longest_common_prefix(first, second, 0, -1, calculated_hashes);
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

    static const size_t CANCELLATION_LENGTH_CACHE_SIZE = 100;

    static Vertex reduce_narrow_slp(
        const Vertex& vertex,
        Cache* calculated_hashes,
        std::unordered_map<Vertex, Vertex>* reduced_vertices
    ) {
      std::vector<std::pair<LongInteger, size_t>> cancellation_length_cache;
      cancellation_length_cache.reserve(CANCELLATION_LENGTH_CACHE_SIZE);

      size_t current_iteration = 0;

      return base_reduce(
          vertex,
          [calculated_hashes, &cancellation_length_cache, &current_iteration](const Vertex& vertex) -> LongInteger {
            ++current_iteration;
            Vertex left = vertex.left_child().negate();
            Vertex right = vertex.right_child();

            if (!left || !right) {
              return LongInteger(0);
            }

            if (crag::slp::get_sub_slp(left, 0, 1) !=
                crag::slp::get_sub_slp(right, 0, 1)) {
              return LongInteger(0);
            }

            size_t previous_reduction_begin = 0;
            size_t previous_reduction_end = cancellation_length_cache.size();

            while (previous_reduction_begin < previous_reduction_end) {
              size_t split = (previous_reduction_begin + previous_reduction_end) / 2;
              const auto& current_length = cancellation_length_cache[split].first;

              if (current_length <= left.length() && current_length <= right.length() &&
                  get_subvertex_hash(left, 0, current_length, calculated_hashes) ==
                  get_subvertex_hash(right, 0, current_length, calculated_hashes)) {
                if (crag::slp::get_sub_slp(left, current_length, current_length + 1) !=
                    crag::slp::get_sub_slp(right, current_length, current_length + 1)) {
                  cancellation_length_cache.at(split).second = current_iteration;
                  return cancellation_length_cache[split].first;
                }
                previous_reduction_begin = split + 1;
              } else {
                previous_reduction_end = split;
              }
            }

            auto new_length_insert_position = cancellation_length_cache.begin() + previous_reduction_begin;

            if (cancellation_length_cache.size() >= CANCELLATION_LENGTH_CACHE_SIZE) {
              auto oldest_entry = std::min_element(cancellation_length_cache.begin(), cancellation_length_cache.end(),
                  [](const std::pair<const LongInteger, size_t>& first, const std::pair<const LongInteger, size_t>& second) {
                    return first.second < second.second;
                  }
              );

              if (new_length_insert_position > oldest_entry) {
                --new_length_insert_position;
                --previous_reduction_begin;
                --previous_reduction_end;
              }

              cancellation_length_cache.erase(oldest_entry);
            }

            assert(new_length_insert_position <= cancellation_length_cache.end());
            assert(new_length_insert_position >= cancellation_length_cache.begin());

            LongInteger cancellation_length = get_longest_common_prefix(
                left,
                right,
                previous_reduction_begin > 0 ? cancellation_length_cache[previous_reduction_begin - 1].first : 0,
                previous_reduction_end < cancellation_length_cache.size() ? cancellation_length_cache[previous_reduction_end].first : -1,
                calculated_hashes
            );

            auto new_length = cancellation_length_cache.emplace(new_length_insert_position, std::move(cancellation_length), current_iteration);

            return new_length->first;
          },
          reduced_vertices
      );
    }

    static Vertex reduce_narrow_slp(const Vertex& vertex) {
      Cache cache;
      std::unordered_map<Vertex, Vertex> reduced_vertices;
      return reduce_narrow_slp(vertex, &cache, &reduced_vertices);
    }


    //! Remove redundant vertices with the same hash value reducing the size of SLP in this way.
    static Vertex remove_duplicates(const Vertex& root, Cache* p_cache, HashRepresentativesCache* p_hash_cache) {
      assert(p_cache);
      assert(p_hash_cache);

      //returning the canonical vertex corresponding to the vertex hash
      auto map_vertex = [&] (const Vertex& v, const std::unordered_map<Vertex, Vertex>& new_vertices) -> Vertex {
        //side effect of this call is the cache filling for all the subvertices
        auto hash = get_vertex_hash(v, p_cache);
        auto item = p_hash_cache->find(hash);
        if (item != p_hash_cache->end()) {
          return item->second;
        } else {
          if (v.height() == 1) {
            p_hash_cache->insert(std::make_pair(hash, v));
            return v;
          } else {

            auto left_hash = get_vertex_hash(v.left_child(), p_cache);
            auto right_hash = get_vertex_hash(v.right_child(), p_cache);

            auto left_item = p_hash_cache->find(left_hash);
            auto right_item = p_hash_cache->find(right_hash);
            assert(left_item != p_hash_cache->end());
            assert(right_item != p_hash_cache->end());

            NonterminalVertex new_v(left_item->second, right_item->second);
            p_hash_cache->insert(std::make_pair(hash, new_v));

            return Vertex(new_v);
          }
        }
      };


      std::unordered_map<Vertex, Vertex> new_vertices;
      //mapping vertices to new ones so there is a unqiue vertice for each appearing hash value
      map_vertices(root, &new_vertices, map_vertex);

      auto root_item = new_vertices.find(root);
      assert(root_item != new_vertices.end());
      return root_item->second;
    }
};

} //namespace slp
} //namespace crag


namespace std {
template <class... Hashers>
struct hash<crag::slp::TVertexHash<Hashers...>> {
  public:
    size_t operator()(const crag::slp::TVertexHash<Hashers...>& hash) const {
      return hash.get_std_hash();
    }
};
}

#endif /* SLP_VERTEX_HASH_H_ */
