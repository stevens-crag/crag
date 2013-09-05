/*
 * EndomorphsimSLP.h
 *
 *  Created on: Jan 29, 2013
 *      Author: pmorar
 */

#ifndef CRAG_FREE_GROUPS_ENDOMORPHISM_SLP_H_
#define CRAG_FREE_GROUPS_ENDOMORPHISM_SLP_H_

#include <map>
#include <unordered_map>
#include <unordered_set>
#include <random>
#include <algorithm>
#include <functional>
#include <assert.h>
#include <chrono>
#include "slp.h"

namespace crag {


namespace endomorphism_default_parameters {
  //! Default hash algorithms used for free reduction.
  typedef crag::slp::TVertexHashAlgorithms<
    crag::slp::hashers::SinglePowerHash,
    crag::slp::hashers::PermutationHash<crag::Permutation16>
  > WeakReducedVertexHashAlgorithms;

  //! Default hash algorithms used for duplicates removal.
  typedef crag::slp::TVertexHashAlgorithms<
    crag::slp::hashers::ImageLengthHash,
    crag::slp::hashers::SinglePowerHash,
    crag::slp::hashers::PermutationHash<crag::Permutation16>
  > WeakVertexHashAlgorithms;
}

template<typename Automorphism>
class AutomorphismDescription;

/**
 * Represents a free group endomorphism using straight-line programs.
 * @tparam TerminalSymbol terminal symbols class
 */
template<typename TerminalSymbol = int>
class EndomorphismSLP {
public:

  typedef typename std::map<TerminalSymbol, slp::Vertex>::size_type size_type;
  typedef typename std::map<TerminalSymbol, slp::Vertex>::iterator iterator;
  typedef typename std::map<TerminalSymbol, slp::Vertex>::const_iterator const_iterator;
  typedef typename std::map<TerminalSymbol, slp::Vertex>::value_type symbol_image_pair_type;


  //! Creates an identity automorphism.
  EndomorphismSLP() {}

  //use default copy/move constructors/assignments

  //! Returns the identity automorphism
  static EndomorphismSLP identity() {
    return EndomorphismSLP();
  }

  //! Returns the automorphism inverting the specified terminal symbbol.
  /**
   * @param symbol inverted terminal symbol, must be > 0
   */
  static EndomorphismSLP inverter(const TerminalSymbol& symbol) {
    assert(is_positive_terminal_symbol(symbol));
    return EndomorphismSLP(symbol);
  }

  //! Returns the automorphism mapping #symbol to the product #symbol #right_multiplier.
  /**
   * @param symbol mapped terminal symbol, must be > 0
   * @param	right_multiplier	right multiplier
   */
  static EndomorphismSLP right_multiplier(const TerminalSymbol& symbol, const TerminalSymbol& right_multiplier) {
    assert(is_positive_terminal_symbol(symbol));
    assert(symbol != right_multiplier);
    EndomorphismSLP tmp;
    tmp.images_.insert(std::make_pair(symbol,
      slp::NonterminalVertex(TerminalVertex(symbol), TerminalVertex(right_multiplier))));
    return tmp;
  }

  //! Returns the automorphism mapping #symbol to the product #left_multiplier #symbol.
  /**
   * @param  left_multiplier  left multiplier
   * @param symbol mapped terminal symbol, must be > 0
   */
  static EndomorphismSLP left_multiplier(const TerminalSymbol& left_multiplier, const TerminalSymbol& symbol) {
    assert(is_positive_terminal_symbol(symbol));
    assert(left_multiplier != symbol);
    EndomorphismSLP tmp;
    tmp.images_.insert(std::make_pair(symbol,
      slp::NonterminalVertex(TerminalVertex(left_multiplier), TerminalVertex(symbol))));
    return tmp;
  }

  //! Returns the automorphism mapping #multiplied to its product with #multiplier. Flag #is_left_multiplication sets the order.
  static EndomorphismSLP multiplication(const TerminalSymbol& multiplied, const TerminalSymbol& multiplier, bool is_left_multiplication) {
    assert(multiplied != multiplier);
    assert(multiplied != -multiplier);
    EndomorphismSLP tmp;
    if (is_left_multiplication) {
      tmp.images_.insert(std::make_pair(multiplied,
        slp::NonterminalVertex(TerminalVertex(multiplier), TerminalVertex(multiplied))));
    } else {
      tmp.images_.insert(std::make_pair(multiplied,
        slp::NonterminalVertex(TerminalVertex(multiplied), TerminalVertex(multiplier))));
    }
    return tmp;
  }

  //! Applies provided function to each inverter, left and right multiplier for the given rank
  /**
   * Enumerates all inverters, left and right miltiplier morphisms for the given rank and apply
   * to them the provided function.
   * @tparam Func type of function applied to morphisms
   */
  template<typename Func>
  static void for_each_basic_morphism(int rank, Func f);

  //! Applies provided function to each multiplication that can be obtained for the given rank.
  template<typename Func>
  static void for_each_multiplication(int rank, Func f);

  //! Returns the composition of endomorphisms specified by the range
  template<typename Iterator>
  static EndomorphismSLP composition(Iterator begin, Iterator end) {
    return identity().compose_with(begin, end);
  }

  //! Returns the composition of #num endomorphisms produced by #generator
  /**
   * @param num      the number of endomorphisms to compose
   * @param generator endomorphisms generator (using operator())
   */
  template<typename Generator>
  static EndomorphismSLP composition(unsigned int num, Generator& generator) {
    return identity().compose_with(num, generator);
  }

  //! Compose with the given endomorphism.
  EndomorphismSLP& operator*=(const EndomorphismSLP& a);

  //! Compose with endomorphisms specified by the range.
  template<typename Iterator>
  EndomorphismSLP& compose_with(Iterator begin, Iterator end) {
    for(;begin != end; ++begin)
      (*this) *= *begin;
    return *this;
  }

  //! Compose with #num automorphism produced by #generator
  /**
   * @param num the number of endomorphisms to compose with
   * @param generator endomorphisms generator (using operator())
   */
  template<typename Generator>
  EndomorphismSLP& compose_with(unsigned int num, Generator& generator) {
    for (unsigned int i = 0; i < num; ++i)
      (*this) *= generator();
    return *this;
  }

  //! Conjugate with the endomorphisms from the given range.
  /**
   * @param num the number of endomorphisms to compose with
   * @param generator endomorphisms generator (using operator())
   * @return conjugator * this * conjugator^{-1}
   */
  template<typename Iterator>
  EndomorphismSLP conjugate_with(Iterator begin, Iterator end) const {
    std::vector<EndomorphismSLP> inverses;
    EndomorphismSLP conjugator;
    for(;begin != end; ++begin) {
      const auto& endomorphism = *begin;
      conjugator *= endomorphism;
      inverses.push_back(endomorphism.inverse());
    }
    return conjugator * (*this) * EndomorphismSLP::composition(inverses.rbegin(), inverses.rend());
  }

  //! Conjugate with the automorphisms generated by the given generator.
  /**
   * @param num the number of endomorphisms to compose with
   * @param generator endomorphisms generator (using operator())
   * @return a this a^{-1}
   */
  template<typename Generator>
  EndomorphismSLP conjugate_with(unsigned int num, Generator& generator) const {
    std::vector<EndomorphismSLP> inverses;
    inverses.reserve(num);
    EndomorphismSLP conjugator;
    for (unsigned int i = 0; i < num; ++i) {
      auto endomorphism = generator();
      conjugator *= endomorphism;
      inverses.push_back(endomorphism.inverse());
    }
    return conjugator * (*this) * EndomorphismSLP::composition(inverses.rbegin(), inverses.rend());
  }

  //! Conjugate with the automorphisms given by the description.
  EndomorphismSLP conjugate_with(const AutomorphismDescription<EndomorphismSLP>& conjugator) const {
    return conjugator() * (*this) * conjugator.inverse();
  }

  //! Returns the automorphisms inverse
  /**
   * Currently supporsts only inverter and left and right multipliers.
   * @return automorphisms inverse
   * @throws std::invalid_argument if can not invert the automorphism
   */
  EndomorphismSLP inverse() const;

  //! Returns the automorphisms with freely reduced images. Might make mistakes but much faster than precise version.
  template<typename VertexHashAlgorithms = endomorphism_default_parameters::WeakReducedVertexHashAlgorithms>
  EndomorphismSLP free_reduction() const {
    typename VertexHashAlgorithms::Cache vertex_hashes;
    EndomorphismSLP result;
    std::unordered_map<slp::Vertex, slp::Vertex> reduced_vertices;

    for_each_non_trivial_image([&result, &vertex_hashes, &reduced_vertices] (const symbol_image_pair_type& pair) {
      auto reduced = VertexHashAlgorithms::reduce_narrow_slp(pair.second, &vertex_hashes, &reduced_vertices);
      //insert if it is not an identity map
      if (reduced.height() != 1 || TerminalVertex(reduced) != pair.first)
        result.images_.insert(std::make_pair(pair.first, reduced));
    });
    return result;
  }

  //! Returns the automorphisms with freely reduced images. It uses matching tables.
  EndomorphismSLP free_reduction_precise() const {
    EndomorphismSLP result;
    slp::MatchingTable mt;
    std::unordered_map<slp::Vertex, slp::Vertex> reduced_vertices;
    for_each_non_trivial_image([&result, &mt, &reduced_vertices] (const symbol_image_pair_type& pair) {
      auto reduced = slp::reduce(pair.second, &mt, &reduced_vertices);
      //insert if it is not an identity map
      if (reduced.height() != 1 || TerminalVertex(reduced) != pair.first)
        result.images_.insert(std::make_pair(pair.first, reduced));
    });
    return result;
  }

  //! Returns the automorphisms, which contains no vertices with the same hash given by template parameter.
  template<typename VertexHashAlgorithms = endomorphism_default_parameters::WeakVertexHashAlgorithms>
  EndomorphismSLP remove_duplicate_vertices() const {
    typename VertexHashAlgorithms::Cache vertex_hashes;
    typename VertexHashAlgorithms::HashRepresentativesCache hash_representatives;
    EndomorphismSLP result;

    for_each_non_trivial_image([&result, &vertex_hashes, &hash_representatives] (const symbol_image_pair_type& pair) {
      auto rd_vertex = VertexHashAlgorithms::remove_duplicates(pair.second, &vertex_hashes, &hash_representatives);
      //insert if it is not an identity map
      if (rd_vertex.height() != 1 || TerminalVertex(rd_vertex) != pair.first)
        result.images_.insert(std::make_pair(pair.first, rd_vertex));
    });
    return result;
  }

  //! Returns the automorphism with its representaiton in normal form.
  EndomorphismSLP normal_form() const;

  //! Returns the image of the terminal.
  slp::VertexWord<TerminalSymbol> image_word(const TerminalSymbol& t) const {
    bool is_positive = is_positive_terminal_symbol(t);
    return slp::VertexWord<TerminalSymbol>(is_positive ? image(t) : image(-t).negate());
  }

  //! Returns the root of the straight-line program representing the positive terminal symbol.
  /**
   * @param t positive terminal symbol
   */
  slp::Vertex image(const TerminalSymbol& t) const {
    assert(is_positive_terminal_symbol(t));
    auto result = images_.find(t);
    if (result == images_.end()) //if it is not in images_, then it is the identity map.
      return TerminalVertex(t);
    else
      return result->second;
  }

  //! Returns the maximal terminal symbol with non-identity image.
  TerminalSymbol max_non_trivial_image_symbol() const {
    if (images_.empty())
      return TerminalSymbol();
    return images_.rbegin()->first;
  }

  //! Returns the number of non-trivial images.
  size_type non_trivial_images_num() const {
    unsigned int n = 0;
    for (auto key_image: images_) {
      TerminalSymbol key = key_image.first;
      const slp::Vertex& image = key_image.second;
      if (image.height() == 1 && TerminalVertex(image).terminal_symbol() == key)
        continue;
      ++n;
    }
    return n;
  }

  //! Returns range (begin, end) of non-trivial images containg pairs (symbol, image)
  std::pair<iterator, iterator> non_trivial_images_range() {
    return std::make_pair(images_.begin(), images_.end());
  }

  //! Returns range (begin, end) of non-trivial images containg pairs (symbol, image)
  std::pair<const_iterator, const_iterator> non_trivial_images_range() const {
    return std::make_pair(images_.begin(), images_.end());
  }

  iterator begin() {
    return images_.begin();
  }

  iterator end() {
    return images_.end();
  }

  const_iterator begin() const {
    return images_.begin();
  }

  const_iterator end() const {
    return images_.end();
  }

  //! Applies the given function to each pair of non-trivial images (symbol, image)
  template<class Function>
  Function for_each_non_trivial_image(Function fn) {
    return std::for_each(images_.begin(), images_.end(), fn);
  }

  //! Applies the given function to each pair of non-trivial images (symbol, image)
  template<class Function>
  Function for_each_non_trivial_image(Function fn) const {
    return std::move(std::for_each(images_.begin(), images_.end(), fn));
  }

  //! Returns true, if the automorphism is identity, and returns false otherwise.
  bool is_identity() const {
    bool is_id = true;
    for_each_non_trivial_image([&is_id](const symbol_image_pair_type& pair) {
            if (is_id)
              is_id = TerminalVertex(pair.first) == slp::reduce(pair.second);
          });
    return is_id;
  }

  bool operator==(const EndomorphismSLP& a) const;

  bool operator!=(const EndomorphismSLP& a) const {
    return !(*this==a);
  }

  void save_to(std::ostream* out) const;
  static EndomorphismSLP load_from(std::istream* in);

  void save_graphviz(std::ostream* out, const std::string& name = "") const;

  //! Prints human readable format
  void print(std::ostream* out) const;

private:
  typedef typename slp::TerminalVertexTemplate<TerminalSymbol> TerminalVertex;

  //! Checks whether the symbol is not inverse.
  static bool is_positive_terminal_symbol(const TerminalSymbol& symbol) {
    static const TerminalSymbol null_symbol = TerminalSymbol();
    return symbol > null_symbol;
  }

  //! Helper method that copies #vertex such that the children of the new vertex are images of the children of #vertex
  /**
   * It assumes that if #vertex is non-terminal, then #images contains its children images
   */
  slp::Vertex map_vertex(const slp::Vertex& vertex, const std::unordered_map<slp::Vertex, slp::Vertex>& images) const;

  EndomorphismSLP(const TerminalSymbol& inverted) {
    images_.insert(std::make_pair(inverted,
      TerminalVertex(inverted).negate()));
  }

  //! The representation of images of terminal symbols by straight-line programs
  /**
   * If there is no entry for a given terminal symbol, then its image is the terminal itself.
   * Also we use {@code std::map} instead of {@code std::unordered_map} because we would like
   * to iterate over the keys in the specific order defined by operator < for TerminalSymbol.
   */
  std::map<TerminalSymbol, slp::Vertex> images_;
};


//---------------------------------------------------------------------
// Helper functions.

//! Compose given endomorphisms.
/**
 * Returns const to behave consistently with built-in types.
 */
template<typename TerminalSymbol>
EndomorphismSLP<TerminalSymbol> operator*(const EndomorphismSLP<TerminalSymbol>& e1, const EndomorphismSLP<TerminalSymbol>& e2) {
  return EndomorphismSLP<TerminalSymbol>(e1) *= e2;
}


//! Find the maximal height of SLPs, representing the endomorphism
template<typename TerminalSymbol>
unsigned int height(const EndomorphismSLP<TerminalSymbol>& e);


//! Find the total number of vertices in SLPs, representing the endomorphism
template<typename TerminalSymbol>
unsigned int slp_vertices_num(const EndomorphismSLP<TerminalSymbol>& e);


//! Find the total number of vertices in SLPs, representing the endomorphism
template<typename TerminalSymbol>
unsigned int slp_unique_images_num(const EndomorphismSLP<TerminalSymbol>& e);


//! Returns the map contatining the lengths of non-trivial images (actually some of the images might be trivial).
template<typename TerminalSymbol>
std::map<TerminalSymbol, LongInteger> images_length(const EndomorphismSLP<TerminalSymbol>& e);


//! Automorphisms generator
/**
 * @tparam TerminalSymbol          terminal symbol representation. We suppose it has a constructor taking index
 * @tparam RandomEngine            engine generating uniformly random non-negative numbers. See std::random library documentation.
 * @tparam TerminalSymbolIndexType terminal symbol index type
 */
template <typename TerminalSymbol = int,
  typename RandomEngine = std::default_random_engine>
class UniformAutomorphismSLPGenerator {
public:
  typedef int index_type;

  //! Constructs a generator of automorphisms of the free group of the given rank.
  /**
   * @param rank free group rank > 0
   */
  UniformAutomorphismSLPGenerator(index_type rank)
    : UniformAutomorphismSLPGenerator(1, rank, std::make_shared<RandomEngine>(), nullptr)
  {}

  //! Constructs a generator of automorphisms of the free group of the given rank.
  /**
   * @param rank free group rank > 0
   * @param seed random engine seed for creation of a new one
   */
  explicit UniformAutomorphismSLPGenerator(index_type rank, typename RandomEngine::result_type seed)
    : UniformAutomorphismSLPGenerator(1, rank, ::std::make_shared<RandomEngine>(seed), nullptr)
  {}

  //! Constructs a generator of automorphisms of the free group of the given rank.
  /**
   * @param rank free group rank > 0
   * @param random_engine random engine
   */
  explicit UniformAutomorphismSLPGenerator(index_type rank, RandomEngine* random_engine)
    : UniformAutomorphismSLPGenerator(1, rank, nullptr, random_engine)
  {}

  //! Constructs a generator of automorphisms of the free group with generators from the interval [min_symbol_index, max_symbol_index].
  /**
   * @param min_symbol_index min index of free group generator > 0
   * @param max_symbol_index min index of free group generator > 0
   */
  UniformAutomorphismSLPGenerator(index_type min_symbol_index, index_type max_symbol_index)
    : UniformAutomorphismSLPGenerator(min_symbol_index, max_symbol_index, std::make_shared<RandomEngine>(), nullptr)
  {}

  //! Constructs a generator of automorphisms of the free group with generators from the interval [min_symbol_index, max_symbol_index].
  /**
   * @param min_symbol_index min index of free group generator > 0
   * @param max_symbol_index min index of free group generator > 0
   * @param seed random engine seed for creation of a new one
   */
  explicit UniformAutomorphismSLPGenerator(index_type min_symbol_index, index_type max_symbol_index, typename RandomEngine::result_type seed)
    : UniformAutomorphismSLPGenerator(min_symbol_index, max_symbol_index, ::std::make_shared<RandomEngine>(seed), nullptr)
  {}

  //! Constructs a generator of automorphisms of the free group with generators from the interval [min_symbol_index, max_symbol_index].
  /**
   * @param min_symbol_index min index of free group generator > 0
   * @param max_symbol_index min index of free group generator > 0
   * @param random_engine random engine
   */
  explicit UniformAutomorphismSLPGenerator(index_type min_symbol_index, index_type max_symbol_index, RandomEngine* random_engine)
    : UniformAutomorphismSLPGenerator(min_symbol_index, max_symbol_index, nullptr, random_engine)
  {}

  //! Set the probability to generate an inverter of terminal symbol. By default it corresponds to uniform distribution.
  void set_inverters_probability(double inverters_probability) {
    assert(inverters_probability >= 0 && inverters_probability <= 1);
    inverters_probability_ = inverters_probability;
  }

  //! Generates a random automorphism.
  EndomorphismSLP<TerminalSymbol> operator()() {
    double p = real_distr_(*random_engine_);
    if (p <= inverters_probability_) {//generate an inverter
      index_type val = inverter_distr_(*random_engine_);
      return EndomorphismSLP<TerminalSymbol>::inverter(TerminalSymbol(MIN_SYMBOL_INDEX + val));
    } else {//generate a multiplier
      index_type val = multiplier_distr_(*random_engine_);
      const bool right_multiplier = (val % 2) == 0;
      val >>= 1;
      const int mapped_symbol_index = MIN_SYMBOL_INDEX + ( val % RANK );
      const TerminalSymbol mapped_symbol(mapped_symbol_index);
      const int multiplier_index = MIN_SYMBOL_INDEX + ( val / RANK );
      const TerminalSymbol multiplier(multiplier_index < mapped_symbol_index ? multiplier_index : multiplier_index + 1);

      if (right_multiplier)
        return EndomorphismSLP<TerminalSymbol>::right_multiplier(mapped_symbol, multiplier);
      else
        return EndomorphismSLP<TerminalSymbol>::left_multiplier(multiplier, mapped_symbol);
    }
  }


private:
  const index_type MIN_SYMBOL_INDEX;
  const index_type MAX_SYMBOL_INDEX;
  const index_type RANK;
  const index_type RIGHT_MULTIPLIERS_COUNT;
  const index_type LEFT_MULTIPLIERS_COUNT = RIGHT_MULTIPLIERS_COUNT;
  const index_type MULTIPLIERS_COUNT = RIGHT_MULTIPLIERS_COUNT + LEFT_MULTIPLIERS_COUNT;
  const index_type INVERTERS_COUNT;
  const index_type COUNT = MULTIPLIERS_COUNT + INVERTERS_COUNT;

  const double DEFAULT_INVERTER_PROBABILITY = static_cast<double>(INVERTERS_COUNT) / COUNT;

  UniformAutomorphismSLPGenerator(index_type min_sym_index, index_type max_sym_index,
                                  const ::std::shared_ptr<RandomEngine>& random_engine_ptr,
                                  RandomEngine* random_engine)
      : MIN_SYMBOL_INDEX(min_sym_index),
        MAX_SYMBOL_INDEX(max_sym_index),
        RANK(max_sym_index - min_sym_index + 1),
        RIGHT_MULTIPLIERS_COUNT(RANK * (RANK - 1)),
        INVERTERS_COUNT(RANK),
        random_engine_ptr_(random_engine_ptr),
        random_engine_(random_engine_ptr ? random_engine_ptr.get() : random_engine),
        inverter_distr_(0, INVERTERS_COUNT - 1),
        multiplier_distr_(0, MULTIPLIERS_COUNT - 1),
        real_distr_(),//interval [0, 1)
        inverters_probability_(DEFAULT_INVERTER_PROBABILITY) {
      assert(min_sym_index > 0);
      assert(max_sym_index > 0);
      assert(max_sym_index >= min_sym_index);
    }

  ::std::shared_ptr<RandomEngine> random_engine_ptr_;
  RandomEngine* random_engine_;
  std::uniform_int_distribution<index_type> inverter_distr_;
  std::uniform_int_distribution<index_type> multiplier_distr_;
  std::uniform_real_distribution<double> real_distr_;
  double inverters_probability_;
};


//! Automorphism description keeping track of the automorphisms and its inverse.
template<typename Automorphism = EndomorphismSLP<int> >
class AutomorphismDescription {
  public:

    AutomorphismDescription()
      : a_(),
        a_inv_(),
        num_(0) {}

    //! Creates the description of a single autmorphism.
    AutomorphismDescription(Automorphism e)
      : a_(std::move(e)),
        a_inv_(std::move(a_.inverse())),
        num_(1u) {}

    //! Generates a random autmorphism.
    template<typename RandomAutomorphismGenerator>
    AutomorphismDescription(unsigned int size, RandomAutomorphismGenerator& random)
      : num_(size) {
      std::vector<Automorphism> parts;
      parts.reserve(size);
      for (unsigned int i = 0; i < size; ++i) {
        parts.push_back(random());
      }
      a_ = Automorphism::composition(parts.begin(), parts.end());
      std::for_each(parts.rbegin(), parts.rend(), [&] (const Automorphism& a) {
        a_inv_ *= a.inverse();
      });
    }

    //! Returns the automorphism itself.
    const Automorphism& operator()() const {
      return a_;
    }

    //! Returns the automorphism inverse.
    const Automorphism& inverse() const {
      return a_inv_;
    }

    //! Returns description of inverse automorphism.
    AutomorphismDescription inverse_description() const {//TODO optimze so it would be without copying
      return AutomorphismDescription(a_inv_, a_, num_);
    }

    //! Returns description of the composition of automorphisms.
    AutomorphismDescription operator*(const AutomorphismDescription& ad) const {
      return AutomorphismDescription(a_ * ad.a_, ad.a_inv_ * a_inv_, num_ + ad.num_);
    }

    //! Modifies itself to be a description of the composition of itself and the given automorphism.
    AutomorphismDescription& operator*=(const AutomorphismDescription& ad) {
      *this = *this * ad;
      return *this;
    }

    //! Conjugate as conjugator * this * conjugator^{-1}
    AutomorphismDescription conjugate_with(const AutomorphismDescription& conjugator) const {
//      Automorphism a = (conjugator() * a_ * conjugator.inverse()).free_reduction();
//      Automorphism a_inv = (conjugator() * a_inv_ * conjugator.inverse()).free_reduction();
      return conjugator * (*this) * conjugator.inverse_description();
    }

    AutomorphismDescription free_reduction() const {
      return AutomorphismDescription(a_.free_reduction(), a_inv_.free_reduction(), num_);
    }

    AutomorphismDescription normal_form() const {
      return AutomorphismDescription(a_.normal_form(), a_inv_.normal_form(), num_);
    }

    AutomorphismDescription remove_duplicate_vertices() const {
      return AutomorphismDescription(a_.remove_duplicate_vertices(), a_inv_.remove_duplicate_vertices(), num_);
    }

    //! Number of composed morphisms consituting the given one.
    unsigned int composed_num() const {
      return num_;
    }

  private:
    Automorphism a_;//todo shared pointer
    Automorphism a_inv_;
    unsigned long num_;//num of composed automorphisms

    AutomorphismDescription(Automorphism a, Automorphism a_inv, unsigned int num)
      : a_(std::move(a)),
        a_inv_(std::move(a_inv)),
        num_(num) {}

};


//! Conjugate the elements of the given vector by the given conjugator.
template<typename T, typename Conjugator>
static std::vector<T> conjugate_all(const std::vector<T>& morphisms, const Conjugator& conjugator) {
  std::vector<T> result;
  result.reserve(morphisms.size());
  for (int i = 0; i < morphisms.size(); ++i) {
    result.push_back(morphisms[i].conjugate_with(conjugator));
  }
  return result;
}


//! Calculates the arguments commutator.
template<typename T>
T make_commutator(const T& first, const T& second) {
  return first * second * first.inverse_description() * second.inverse_description();
}


//! Reduces endomorphisms of their descriptions with optional logging.
class AutomorphismReducer { //todo replace logging to std::cout by logging facility
    class SizePrinter {
      public:
        template<typename AutDescription>
        void operator()(const AutDescription& aut) const {
          std::cout << slp_vertices_num(aut()) << "," << slp_vertices_num(aut.inverse());
        }

        template<typename TerminalSymbol>
        void operator()(const EndomorphismSLP<TerminalSymbol>& aut) const {
          std::cout << slp_vertices_num(aut);
        }
    };

  public:
    template<typename Aut>
    static Aut reduce(const Aut& aut, bool is_logging) {
      return is_logging ? reduce(aut, SizePrinter()) : reduce(aut);
    }

    template<typename Aut, typename AutLogger>
    static Aut reduce(const Aut& aut, AutLogger aut_logger) {
        std::ostream& out = std::cout;
        out << "|";
        aut_logger(aut);
        out << "|";

        out << " fr ";
        auto fr = perform_action_with_logging(aut_logger, [&aut]() {return aut.free_reduction();});

        out << " rd ";
        auto rd = perform_action_with_logging(aut_logger, [&fr]() {return fr.remove_duplicate_vertices();});

        out << " nf ";
        auto result = perform_action_with_logging(aut_logger, [&rd]() {return rd.normal_form();});

        out << std::endl;
        return result;
    }

    template<typename Aut>
    static Aut reduce(const Aut& aut) {
      auto fr = aut.free_reduction();
      auto rd = fr.remove_duplicate_vertices();;
      return rd.normal_form();
    }

    template<typename Func, typename AutLogger>
    static auto perform_action_with_logging(AutLogger aut_logger, Func f) -> decltype(f()) {
     auto start_time = std::chrono::high_resolution_clock::now();
     auto result = f();
     auto duration_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
           std::chrono::high_resolution_clock::now() - start_time
           );
     std::cout << "|";
     aut_logger(result);
     std::cout << "| " << duration_in_ms.count() << "ms";
     return result;
    }
};

//! Set of commutators for a given pair of automorphisms. We use it only because we can not compute inverses efficiently.
template<typename AutDescription>
class CommutatorSet {
  public:
    CommutatorSet(const AutDescription& first, const AutDescription& second) {
      auto inv_first = first.inverse_description();
      auto inv_second = second.inverse_description();
      comm_.reserve(4);
      comm_.push_back(make_commutator(first, second));
      comm_.push_back(make_commutator(inv_first, second));
      comm_.push_back(make_commutator(first, inv_second));
      comm_.push_back(make_commutator(inv_first, inv_second));
    }

    const AutDescription& get(bool first_inversed, bool second_inversed) const {
      return comm_[(first_inversed ? 1 : 0) + (second_inversed ? 2 : 0)];
    }

    CommutatorSet conjugate_with(const AutDescription& conjugator) const {
      return CommutatorSet(std::move(conjugate_all(comm_, conjugator)));
    }

    CommutatorSet reduce() const {
      std::vector<AutDescription> new_comm;
      new_comm.reserve(comm_.size());
      std::transform(comm_.cbegin(), comm_.cend(),
                     std::back_inserter(new_comm),
                     [] (const AutDescription& ad) {return AutomorphismReducer::reduce(ad);});
      return CommutatorSet(std::move(new_comm));
    }

  private:
    CommutatorSet(std::vector<AutDescription>&& v)
      : comm_(v) {}

    std::vector<AutDescription> comm_;

    CommutatorSet() {}
};

//-------------------------------------------------------------------------------------
// Implementation of EndomorphismSLP methods.


template<typename TerminalSymbol>
template<typename Func>
void EndomorphismSLP<TerminalSymbol>::for_each_basic_morphism(int rank, Func f) {
  assert(rank > 0);
  for (int i = 1; i <= rank; ++i)
    f(EndomorphismSLP<TerminalSymbol>::inverter(i));
  for (int i = 1; i <= rank; ++i)
    for (int j = -static_cast<int>(rank); j <= rank; ++j)
      if (j != i && j != -i && j != 0) {
        f(EndomorphismSLP<TerminalSymbol>::left_multiplier(j, i));
        f(EndomorphismSLP<TerminalSymbol>::right_multiplier(i, j));
      }
}

template<typename TerminalSymbol>
template<typename Func>
void EndomorphismSLP<TerminalSymbol>::for_each_multiplication(int rank, Func f) {
  assert(rank > 0);
  const int lower_bound = -static_cast<int>(rank);
  for (int i = lower_bound; i <= rank; ++i) {
    if (i == 0) {
      continue;
    }
    for (int j = lower_bound; j <= rank; ++j) {
      if (j == 0 || j == i || j == -i) {
        continue;
      }
      f(EndomorphismSLP<TerminalSymbol>::multiplication(i, j, true));
      f(EndomorphismSLP<TerminalSymbol>::multiplication(i, j, false));
    }
  }
}



template <typename TerminalSymbol>
EndomorphismSLP<TerminalSymbol> EndomorphismSLP<TerminalSymbol>::inverse() const {
  if (images_.size() == 0) //identity
    return *this;
  if (images_.size() > 1)
    throw std::invalid_argument("Unsupported endomorphism with more than one non-trivial terminal image!");

  auto img = images_.begin();
  TerminalSymbol symbol = img->first;
  slp::Vertex image = img->second;
  if (image.height() > 2)
    throw std::invalid_argument("Unsupported endomorphism with unsupported slp height > 2!");

  if (image.height() == 1) {//inverter
    return *this;
  }

  TerminalVertex left(image.left_child());
  TerminalVertex right(image.right_child());

  TerminalSymbol left_symbol = left.terminal_symbol();
  TerminalSymbol right_symbol = right.terminal_symbol();

  if (! (left_symbol == symbol || right_symbol == symbol))
    throw std::invalid_argument("Unsupported endomorphism not mapping the symbol to the product of another one and itself!");

  if (left_symbol == symbol) {
    return right_multiplier(symbol, -right_symbol);
  } else {
    return left_multiplier(-left_symbol, symbol);
  }
}

template <typename TerminalSymbol>
bool EndomorphismSLP<TerminalSymbol>::operator==(const EndomorphismSLP<TerminalSymbol>& a) const {
  if (this == &a)
    return true;
  if (non_trivial_images_num() != a.non_trivial_images_num())
    return false;
  auto images = non_trivial_images_range();
  auto a_images = a.non_trivial_images_range();

  //checking that the same letters have non-trivial images
  auto img_iterator = images_.begin();
  auto a_img_iterator = a.images_.begin();
  while (img_iterator != images_.end()) {
    if (img_iterator->first != a_img_iterator->first)
      return false;
    ++img_iterator;
    ++a_img_iterator;
  }

  //checking that the images are the same
  slp::MatchingTable mt;
  img_iterator = images_.begin();
  a_img_iterator = a.images_.begin();
  while (img_iterator != images_.end()) {
    slp::VertexWord<TerminalSymbol> word(img_iterator->second);
    slp::VertexWord<TerminalSymbol> a_word(a_img_iterator->second);
    if (!word.is_equal_to(a_word, &mt))
      return false;
    ++img_iterator;
    ++a_img_iterator;
  }

  return true;
}

template <typename TerminalSymbol>
EndomorphismSLP<TerminalSymbol>& EndomorphismSLP<TerminalSymbol>::operator*=(const EndomorphismSLP<TerminalSymbol>& a) {
  std::unordered_map<slp::Vertex, slp::Vertex> new_vertices;//a's vertices to new vertices correspondence

  for (const auto& root_entry: a.images_) {//mapping vertices of #a to new ones
    slp::map_vertices(root_entry.second, &new_vertices,
                      std::bind(&EndomorphismSLP<TerminalSymbol>::map_vertex, *this, std::placeholders::_1, std::placeholders::_2));
  }

  //replacing roots
  std::map<TerminalSymbol, slp::Vertex> new_images;
  for (const auto& root_entry: a.images_) {
    auto new_root = new_vertices.find(root_entry.second)->second;
    new_images.insert(std::make_pair(root_entry.first, new_root));
  }
  //adding images that were not inserted
  for (const auto& root_entry: images_) {
    if (new_images.find(root_entry.first) == new_images.end())//it was not mapped by a
      new_images.insert(root_entry);
  }

  swap(images_, new_images);
  return *this;
}

template<typename TerminalSymbol>
slp::Vertex EndomorphismSLP<TerminalSymbol>::map_vertex(const slp::Vertex& vertex, const std::unordered_map<slp::Vertex, slp::Vertex>& images) const {
  auto item = images.find(vertex.negate());
  if (item != images.end()) {//already mapped inverse
    return item->second.negate();
  }

  if (!vertex)
    return vertex;//Mapping null vertex

  if (vertex.height() == 1) {//the vertex is terminal so map it to our corresponding root
    TerminalSymbol symbol = TerminalVertex(vertex).terminal_symbol();
    bool is_positive = is_positive_terminal_symbol(symbol);
    TerminalSymbol positive_symbol = is_positive ? symbol : -symbol;
    slp::Vertex v = image(positive_symbol);
    if (TerminalVertex(positive_symbol) == v)//if id map
      return vertex;
    return is_positive ? v : v.negate();
  } else {//for a nonterminal we already processed its children because postorder inspector
    auto left_val = images.find(vertex.left_child());
    auto right_val = images.find(vertex.right_child());
    assert(left_val != images.end() && right_val != images.end());
    const slp::Vertex& left = left_val->second;
    const slp::Vertex& right = right_val->second;
    if (left == vertex.left_child() && right == vertex.right_child()) //if children were not copied, then we should not copy vertex
      return vertex;
    return slp::NonterminalVertex(left, right);
  }
}

template<typename TerminalSymbol>
void EndomorphismSLP<TerminalSymbol>::save_to(std::ostream* out) const {
  long vertex_num = 0;

  std::unordered_map<size_t, std::pair<long, long>> non_terminals;
  std::unordered_map<size_t, TerminalSymbol> terminals;

  //we save the order in wich vertices occur during postorder inspection because we want to save them in this
  //order so it would be easier to reconstruct the automorphism later
  std::vector<long> non_terminals_order;
  std::vector<long> terminals_order;

  auto processor = [&] (const slp::Vertex& vertex, std::unordered_map<slp::Vertex, long>& mapped_images) {
    auto item = mapped_images.find(vertex.negate());
    if (item != mapped_images.end()) {//inverse was already visited
      return -item->second; //negate_index
    }

    ++vertex_num;
    if (vertex.height() == 1) {//the vertex is terminal
      TerminalSymbol symbol = TerminalVertex(vertex).terminal_symbol();
      bool is_positive = symbol > 0;
      const TerminalSymbol positive_symbol = is_positive ? symbol : - symbol;
      terminals.insert(std::make_pair(vertex_num, positive_symbol));
      terminals_order.push_back(vertex_num);
      if (is_positive) {
        mapped_images.insert(std::make_pair(vertex.negate(), -vertex_num));
        return vertex_num;
      } else {
        mapped_images.insert(std::make_pair(vertex.negate(), vertex_num));
        return -vertex_num;
      }
    } else {//nonterminal
      size_t left_val = mapped_images.find(vertex.left_child())->second;
      size_t right_val = mapped_images.find(vertex.right_child())->second;
      non_terminals.insert(std::make_pair(vertex_num, std::make_pair(left_val, right_val)));
      non_terminals_order.push_back(vertex_num);
    }
    return vertex_num;
  };

  std::unordered_map<slp::Vertex, long> vertex_numbers;
  for (const auto& root_entry: images_) {
    slp::map_vertices(root_entry.second, &vertex_numbers,
                      processor);
  }

  //writing
  //"number of nontrivial images" "number of terminals" "number of non-terminals"
  *out << images_.size() << " " << terminals.size() << " " << non_terminals.size() << std::endl;

  //writing terminal symbols
  //"terminal vertex index" "terminal symbol"
  for (long terminal_vertex_index: terminals_order) {
    const auto& terminal = *(terminals.find(terminal_vertex_index));
    *out << terminal.first << " " << terminal.second << std::endl;
  }

  //writing non-terminals
  //"vertex index" "left child index" "right child index"
  for (size_t non_terminal_vertex_index: non_terminals_order) {
    const auto& non_terminal = *(non_terminals.find(non_terminal_vertex_index));
    *out << non_terminal.first << " " << non_terminal.second.first << " " << non_terminal.second.second << std::endl;
  }

  //writing roots
  //"terminal symbol" "vertex index"
  for (const auto& root_entry: images_)
    *out << root_entry.first << " " << vertex_numbers.find(root_entry.second)->second << std::endl;
}

template<typename TerminalSymbol>
EndomorphismSLP<TerminalSymbol> EndomorphismSLP<TerminalSymbol>::load_from(std::istream* in) {
  size_t roots_num;
  size_t terminals_num;
  size_t non_terminals_num;
  *in >> roots_num >> terminals_num >> non_terminals_num;

  std::unordered_map<long, slp::Vertex> vertices;
  for (int i = 0; i < terminals_num; ++i) {
    long index;
    TerminalSymbol image;
    *in >> index >> image;
    in->ignore();
    vertices.insert(std::make_pair(index, TerminalVertex(image)));
  }

  auto get_vertex = [&vertices] (long index) {
    bool is_positive = index > 0;
    long positive_index = is_positive ? index : -index;
    const slp::Vertex& v = vertices.find(positive_index)->second;
    return is_positive ? v : v.negate();
  };

  for (int i = 0; i < non_terminals_num; ++i) {
    long index;
    long l_index;
    long r_index;
    *in >> index >> l_index >> r_index;
    in->ignore();
    const slp::Vertex& left = get_vertex(l_index);
    const slp::Vertex& right = get_vertex(r_index);
    vertices.insert(std::make_pair(index, slp::NonterminalVertex(left, right)));
  }

  EndomorphismSLP e;
  for (int i = 0; i < roots_num; ++i) {
    TerminalSymbol key;
    size_t index;
    *in >> key >> index;
    in->ignore();
    const slp::Vertex& root = get_vertex(index);
    e.images_.insert(std::make_pair(key, root));
  }
  return e;
}


template<typename TerminalSymbol>
void EndomorphismSLP<TerminalSymbol>::save_graphviz(std::ostream *p_out, const std::string& name) const {
  static const char* INDENT = "\t";

  std::ostream& out = *p_out;
  out << "digraph " << name << " {" << std::endl;
  out << "node [shape=point]" << std::endl;

  long vertex_num = 0;

  std::unordered_map<size_t, std::pair<long, long>> non_terminals;
  std::unordered_map<size_t, TerminalSymbol> terminals;
  std::unordered_map<TerminalSymbol, size_t> sym_to_vertex_num;

  auto processor = [&] (const slp::Vertex& vertex, std::unordered_map<slp::Vertex, long>& mapped_images) {
    auto item = mapped_images.find(vertex.negate());
    if (item != mapped_images.end()) {//inverse was already visited
      auto negate_index = item->second;
      return -negate_index;
    }

    ++vertex_num;
    if (vertex.height() == 1) {//the vertex is terminal
      TerminalSymbol symbol = TerminalVertex(vertex).terminal_symbol();
      bool is_positive = symbol > 0;
      const TerminalSymbol positive_symbol = is_positive ? symbol : - symbol;
      terminals.insert(std::make_pair(vertex_num, positive_symbol));
      sym_to_vertex_num.insert(std::make_pair(positive_symbol, vertex_num));
      if (is_positive) {
        mapped_images.insert(std::make_pair(vertex.negate(), -vertex_num));
        return vertex_num;
      } else {
        mapped_images.insert(std::make_pair(vertex.negate(), vertex_num));
        return -vertex_num;
      }
    } else {//nonterminal
      size_t left_val = mapped_images.find(vertex.left_child())->second;
      size_t right_val = mapped_images.find(vertex.right_child())->second;
      non_terminals.insert(std::make_pair(vertex_num, std::make_pair(left_val, right_val)));
    }
    return vertex_num;
  };

  std::unordered_map<slp::Vertex, long> vertex_numbers;
  for (const auto& root_entry: images_) {
    slp::map_vertices(root_entry.second, &vertex_numbers,
                      processor);
  }

  //writing root styles
  for (const auto& root_entry: images_) {
    out << INDENT << "\"i" << root_entry.first << "\" [shape=plaintext, label=\"" << root_entry.first << "\"];" << std::endl;

    out << INDENT << "\"i" << root_entry.first << "\" -> ";
    const slp::Vertex& img = root_entry.second;
    if (img.height() <= 1) {
      TerminalSymbol symbol = TerminalVertex(img).terminal_symbol();
      bool is_positive = symbol > 0;
      TerminalSymbol positive_symbol = is_positive ? symbol : -symbol;
      const auto& terminal = sym_to_vertex_num.find(positive_symbol);
      if (terminal != sym_to_vertex_num.end()) {
        if (is_positive) {
          out << terminal->second << ";" << std::endl;
        } else {
          out << terminal->second << "[label=\"-\"];" << std::endl;
        }
      } else {
        out << "\"" << symbol << "\"[shape=plaintext];" << std::endl;
      }
    } else {
      auto non_terminal_index = vertex_numbers.find(img)->second;
      if (non_terminal_index > 0) {
        out << non_terminal_index << " [style=dotted];" << std::endl;
      } else {
        out << (-non_terminal_index) << " [style=dotted, label=\"-\"];" << std::endl;
      }
    }
  }

  //writing terminals style
  for (const auto& pair: terminals) {
    out << INDENT << pair.first << " [shape=plaintext, label=\"" << pair.second << "\"];" << std::endl;
  }

  //writing non-terminals
  for (const auto& non_terminal: non_terminals) {
    size_t non_terminal_index = non_terminal.first;

    auto left_index = non_terminal.second.first;
    auto right_index = non_terminal.second.second;
    if (left_index > 0) {
      out << INDENT << non_terminal_index << " -> " << left_index << ";" << std::endl;
    } else {
      out << INDENT << non_terminal_index << " -> " << -left_index << "[label=\"-\"];" << std::endl;
    }
    if (right_index > 0){
      out << INDENT << non_terminal_index << " -> " << right_index << "[color=red,style=dashed];" << std::endl;
    } else {
      out << INDENT << non_terminal_index << " -> " << -right_index << "[color=red,style=dashed,label=\"-\"];" << std::endl;
    }
  }


  out << "}" << std::endl;
}

template<typename TerminalSymbol>
void EndomorphismSLP<TerminalSymbol>::print(std::ostream *p_out) const {
  std::ostream& out = *p_out;

  size_t vertex_num = 0;

  std::unordered_map<size_t, std::pair<size_t, size_t>> non_terminals;
  std::unordered_map<size_t, TerminalSymbol> terminals;

  auto processor = [&] (const slp::Vertex& vertex, const std::unordered_map<slp::Vertex, size_t>& mapped_images) {
    ++vertex_num;
    if (vertex.height() == 1) {//the vertex is terminal
      const TerminalSymbol& symbol = TerminalVertex(vertex).terminal_symbol();
      terminals.insert(std::make_pair(vertex_num, symbol));
    } else {//nonterminal
      size_t left_val = mapped_images.find(vertex.left_child())->second;
      size_t right_val = mapped_images.find(vertex.right_child())->second;
      non_terminals.insert(std::make_pair(vertex_num, std::make_pair(left_val, right_val)));
    }
    return vertex_num;
  };

  std::unordered_map<slp::Vertex, size_t> vertex_numbers;
  for (auto root_entry: images_) {
    slp::map_vertices(root_entry.second, &vertex_numbers,
                      processor);
  }

  std::unordered_map<size_t, TerminalSymbol> non_terminals_to_roots;

  //writing root styles
  for (auto root_entry: images_) {
    auto img = root_entry.second;
    if (img.height() <= 1) {
      auto symbol = TerminalVertex(img).terminal_symbol();
      out << "\"" << root_entry.first << "\" -> \"" << symbol << "\"" << std::endl;
    } else {
      non_terminals_to_roots.insert(std::make_pair(vertex_numbers.find(img)->second, root_entry.first));
    }
  }

  //writing non-terminals
  for (auto non_terminal: non_terminals) {
    size_t non_terminal_index = non_terminal.first;

    size_t left_index = non_terminal.second.first;
    size_t right_index = non_terminal.second.second;
    auto left_terminal = terminals.find(left_index);
    auto right_terminal = terminals.find(right_index);

    auto root_symbol = non_terminals_to_roots.find(non_terminal_index);
    if (root_symbol == non_terminals_to_roots.end()) {
      out << non_terminal_index;
    } else {
      out << "\"" << root_symbol->second << "\"";
    }

    out << " -> ";
    if (left_terminal == terminals.end()) {
      out << left_index;
    } else {
      out << "\"" << left_terminal->second << "\"";
    }
    out << " ";
    if (right_terminal == terminals.end()) {
      out << right_index;
    } else {
      out << "\"" << right_terminal->second << "\"";
    }
    out << std::endl;
  }
}


namespace internal {
  class Packer {//todo remove class
    public:
      static slp::Vertex pack_slps_into_one(const std::vector<slp::Vertex>& v) {
        std::size_t num = v.size();
        assert (num > 0);
        if (num == 1) {
          return v[0];
        }
        std::vector<slp::Vertex> packed_v;
        unsigned int packed_num = (num / 2) + (num % 2);
        packed_v.reserve(packed_num);
        for (std::size_t i = 0; i < num - 1; i += 2) {
          packed_v.push_back(slp::NonterminalVertex(v[i], v[i + 1]));
        }
        if (num % 2 == 1) {
          packed_v.push_back(v.back());
        }
        return pack_slps_into_one(packed_v);
      }
  };
}

template <typename TerminalSymbol>
EndomorphismSLP<TerminalSymbol> EndomorphismSLP<TerminalSymbol>::normal_form() const {
  if (non_trivial_images_num() == 0)
    return EndomorphismSLP();
  //we rewrite all vertices into a single SLP then find normal form
  //and then split into pieces according to original vertices lengths

  //end positions in the unifying SLP
  std::vector<std::pair<TerminalSymbol, LongInteger> > end_positions;

  //vertices to combine into SLP
  std::vector<slp::Vertex> vertices;
  vertices.reserve(non_trivial_images_num());

  LongInteger length_accumulator(0);
  for_each_non_trivial_image([&] (const symbol_image_pair_type& pair) {
    const TerminalSymbol& k = pair.first;
    const slp::Vertex& v = pair.second;

    //fill positions
    length_accumulator += v.length();
    end_positions.push_back(std::make_pair(k, length_accumulator));

    //adding vertices to process
    vertices.push_back(v);
  });

  slp::Vertex slp = internal::Packer::pack_slps_into_one(vertices);
  auto nf = slp::recompression::normal_form(slp);

  //extracting slps
  EndomorphismSLP result;
  LongInteger begin(0);
  for (auto key_pos_pair: end_positions) {
    auto end = key_pos_pair.second;
    auto sub_slp = slp::get_sub_slp(nf, begin, end);
    result.images_.insert(std::make_pair(key_pos_pair.first, sub_slp));
    begin = end;
  }

  return result;
}

//-------------------------------------------------------------------------
// Implementation of helper functions.


template<typename TerminalSymbol>
unsigned int height(const EndomorphismSLP<TerminalSymbol>& e) {
  unsigned int h = 0;
  auto pick_max_height = [&h] (const typename EndomorphismSLP<TerminalSymbol>::symbol_image_pair_type& v) {
    const unsigned int v_h = v.second.height();
    if (v_h > h)
      h = v_h;
  };
  e.for_each_non_trivial_image(pick_max_height);
  return h;
}


template<typename TerminalSymbol>
unsigned int slp_vertices_num(const EndomorphismSLP<TerminalSymbol>& e) {
  std::unordered_set<slp::Vertex> visited_vertices;

  auto acceptor = [&visited_vertices] (const slp::inspector::InspectorTask& task) {
    return visited_vertices.count(task.vertex) == 0
        && visited_vertices.count(task.vertex.negate()) == 0;
  };

  auto inspect_root =[&acceptor,&visited_vertices] (const typename EndomorphismSLP<TerminalSymbol>::symbol_image_pair_type& v) {
    slp::Inspector<slp::inspector::Postorder, decltype(acceptor)> inspector(v.second, acceptor);
    while (!inspector.stopped()) {
      visited_vertices.insert(inspector.vertex());
      inspector.next();
    }
  };

  e.for_each_non_trivial_image(inspect_root);
  return visited_vertices.size();
}

template<typename TerminalSymbol>
unsigned int slp_unique_images_length_num(const EndomorphismSLP<TerminalSymbol>& e) {
  std::set<LongInteger> visited_vertices;//map sum images length -> vertex

  auto acceptor = [&visited_vertices] (const slp::inspector::InspectorTask& task) {
    return visited_vertices.find(task.vertex.length()) == visited_vertices.end();
  };

  auto inspect_root =[&acceptor,&visited_vertices] (const typename EndomorphismSLP<TerminalSymbol>::symbol_image_pair_type& v) {
    slp::Inspector<slp::inspector::Postorder, decltype(acceptor)> inspector(v.second, acceptor);
    while (!inspector.stopped()) {
      visited_vertices.insert(inspector.vertex().length());
      inspector.next();
    }
  };

  e.for_each_non_trivial_image(inspect_root);
  return visited_vertices.size();
}


template<typename TerminalSymbol>
std::map<TerminalSymbol, LongInteger> images_length(const EndomorphismSLP<TerminalSymbol>& e) {
  std::map<TerminalSymbol, LongInteger> key_to_lengths;
  auto add_length = [&key_to_lengths] (const typename EndomorphismSLP<TerminalSymbol>::symbol_image_pair_type& pair) {
   auto key = pair.first;
   slp::Vertex v = pair.second;
   key_to_lengths.insert(std::make_pair(key, v.length()));
  };
  e.for_each_non_trivial_image(add_length);
  return key_to_lengths;
}

} /* namespace crag */
#endif /* CRAG_FREE_GROUPS_ENDOMORPHISM_SLP_H_ */
