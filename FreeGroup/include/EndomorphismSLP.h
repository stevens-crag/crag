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
 */
class EndomorphismSLP {
public:
  typedef slp::TerminalSymbol TerminalSymbol;
  typedef slp::TerminalVertex TerminalVertex;

  typedef std::map<TerminalSymbol, slp::Vertex>::size_type size_type;
  typedef std::map<TerminalSymbol, slp::Vertex>::iterator iterator;
  typedef std::map<TerminalSymbol, slp::Vertex>::const_iterator const_iterator;
  typedef std::map<TerminalSymbol, slp::Vertex>::value_type symbol_image_pair_type;


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

  //! Compose with the given endomorphism.
  EndomorphismSLP operator*(const EndomorphismSLP& a) const {
    return EndomorphismSLP(*this) *= a;
  }


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
  EndomorphismSLP conjugate_with(const AutomorphismDescription<EndomorphismSLP>& conjugator) const;
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
    auto reducer = [&vertex_hashes] (const slp::Vertex& vertex,
        std::unordered_map<slp::Vertex, slp::Vertex>* p_reduced_vertices) {
      return VertexHashAlgorithms::reduce_narrow_slp(vertex, &vertex_hashes, p_reduced_vertices);
    };

    return free_reduction_internal(&reducer);
  }

  //! Returns the automorphisms with freely reduced images. It uses matching tables.
  EndomorphismSLP free_reduction_precise() const {
    slp::MatchingTable mt;
    auto reducer = [&mt] (const slp::Vertex& vertex,
        std::unordered_map<slp::Vertex, slp::Vertex>* p_reduced_vertices) {
      return slp::reduce(vertex, &mt, p_reduced_vertices);
    };

    return free_reduction_internal(&reducer);
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
  slp::VertexWord image_word(const TerminalSymbol& t) const {
    bool is_positive = is_positive_terminal_symbol(t);
    return slp::VertexWord(is_positive ? image(t) : image(-t).negate());
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

  template<typename Reducer>
  EndomorphismSLP free_reduction_internal(Reducer* p_reducer) const {
    EndomorphismSLP result;
    std::unordered_map<slp::Vertex, slp::Vertex> reduced_vertices;
    for_each_non_trivial_image([&result, &reduced_vertices, p_reducer] (const symbol_image_pair_type& pair) {
      auto reduced = p_reducer->operator()(pair.second, &reduced_vertices);
      //insert if it is not an identity map
      if (reduced.height() != 1 || TerminalVertex(reduced) != pair.first)
        result.images_.insert(std::make_pair(pair.first, reduced));
    });
    return result;
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

//! Find the maximal height of SLPs, representing the endomorphism
unsigned int height(const EndomorphismSLP& e);


//! Find the total number of vertices in SLPs, representing the endomorphism
unsigned int slp_vertices_num(const EndomorphismSLP& e);


//! Find the total number of vertices in SLPs, representing the endomorphism
unsigned int slp_unique_images_num(const EndomorphismSLP& e);


//! Returns the map contatining the lengths of non-trivial images (actually some of the images might be trivial).
std::map<slp::TerminalSymbol, LongInteger> images_length(const EndomorphismSLP& e);


//! Automorphisms generator
/**
 * @tparam RandomEngine            engine generating uniformly random non-negative numbers. See std::random library documentation.
 */
template <
  typename RandomEngine = std::default_random_engine
>
class UniformAutomorphismSLPGenerator {
public:
  typedef size_t index_type;
  typedef slp::TerminalSymbol TerminalSymbol;

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
  UniformAutomorphismSLPGenerator(index_type rank, typename RandomEngine::result_type seed)
    : UniformAutomorphismSLPGenerator(1, rank, ::std::make_shared<RandomEngine>(seed), nullptr)
  {}

  //! Constructs a generator of automorphisms of the free group of the given rank.
  /**
   * @param rank free group rank > 0
   * @param random_engine random engine
   */
  UniformAutomorphismSLPGenerator(index_type rank, RandomEngine* random_engine)
    : UniformAutomorphismSLPGenerator(1, rank, nullptr, random_engine)
  {}

  //! Constructs a generator of automorphisms of the free group with generators from the interval [min_symbol_index, max_symbol_index].
  /**
   * @param min_symbol_index min index of free group generator > 0
   * @param max_symbol_index min index of free group generator > 0
   */
  //UniformAutomorphismSLPGenerator(index_type min_symbol_index, index_type max_symbol_index)
  //  : UniformAutomorphismSLPGenerator(min_symbol_index, max_symbol_index, std::make_shared<RandomEngine>(), nullptr)
  //{}

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
  EndomorphismSLP operator()() {
    double p = real_distr_(*random_engine_);
    if (p <= inverters_probability_) {//generate an inverter
      index_type val = inverter_distr_(*random_engine_);
      return EndomorphismSLP::inverter(TerminalSymbol(MIN_SYMBOL_INDEX + val));
    } else {//generate a multiplier
      index_type val = multiplier_distr_(*random_engine_);
      const bool right_multiplier = (val % 2) == 0;
      val >>= 1;
      const int mapped_symbol_index = MIN_SYMBOL_INDEX + ( val % RANK );
      const TerminalSymbol mapped_symbol(mapped_symbol_index);
      const int multiplier_index = MIN_SYMBOL_INDEX + ( val / RANK );
      const TerminalSymbol multiplier(multiplier_index < mapped_symbol_index ? multiplier_index : multiplier_index + 1);

      if (right_multiplier)
        return EndomorphismSLP::right_multiplier(mapped_symbol, multiplier);
      else
        return EndomorphismSLP::left_multiplier(multiplier, mapped_symbol);
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
template<typename Automorphism = EndomorphismSLP>
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

    //! Returns the automorphism itself.
    const Automorphism& aut() const {
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

        void operator()(const EndomorphismSLP& aut) const {
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
// Implementation of templated EndomorphismSLP methods.


template<typename Func>
void EndomorphismSLP::for_each_basic_morphism(int rank, Func f) {
  assert(rank > 0);
  for (int i = 1; i <= rank; ++i)
    f(EndomorphismSLP::inverter(i));
  for (int i = 1; i <= rank; ++i)
    for (int j = -static_cast<int>(rank); j <= rank; ++j)
      if (j != i && j != -i && j != 0) {
        f(EndomorphismSLP::left_multiplier(j, i));
        f(EndomorphismSLP::right_multiplier(i, j));
      }
}

template<typename Func>
void EndomorphismSLP::for_each_multiplication(int rank, Func f) {
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
      f(EndomorphismSLP::multiplication(i, j, true));
      f(EndomorphismSLP::multiplication(i, j, false));
    }
  }
}





} /* namespace crag */
#endif /* CRAG_FREE_GROUPS_ENDOMORPHISM_SLP_H_ */
