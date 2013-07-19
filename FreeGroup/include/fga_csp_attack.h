#ifndef CRAG_FREE_GROUPS_FGA_CSP_ATTACK_H
#define CRAG_FREE_GROUPS_FGA_CSP_ATTACK_H

/*
 * fga_csp_attack.h
 *
 *  Created on: July 19, 2013
 *      Author: pmorar
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <functional>


#include "EndomorphismSLP.h"

namespace crag {
namespace fga_csp_attack {


typedef std::chrono::high_resolution_clock our_clock;

//----------------------------------------------------------------------------------------
// Statistic accumulators

//! Our statistics accumulator.
template<typename T>
class Statistic {
  public:
    void operator()(const T& val) {
      if (count_++ > 0) {
        sum_ += val;
        min_ = min_ > val ? val : min_;
        max_ = max_ < val ? val : max_;
      } else {
        min_ = max_ = sum_ = val;
      }
    }

    T max() const {
      return max_;
    }

    T min() const {
      return min_;
    }

    T sum() const {
      return sum_;
    }

    T count() const {
      return count_;
    }

  private:
    T max_;
    T min_;
    T sum_;
    T count_;
};

template<typename Stat>
auto average(const Stat& stat) -> decltype(stat.sum() / stat.count()) {
  return stat.sum() / stat.count();
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const Statistic<T>& stat) {
  return out << "(avg = " << average(stat) << ", "
                 << "min = " << stat.min() << ", "
                 << "max = " << stat.max() << ")";
}

//----------------------------------------------------
// Helper methods

//using Endomorphism with alphabet represented by integers.
typedef EndomorphismSLP<int> Aut;

//! Accumulates statistics for images lengths
Statistic<LongInteger> get_endomorphism_images_lengths_stat(const Aut& e);

//! Returns the sum of lengths of generators images
LongInteger total_images_length(const Aut& em);

//! Returns the max of lengths of generators images
LongInteger max_images_length(const Aut& em);

//! Returns the square root of sum of squares of images length
LongInteger sqrt_of_sq_of_img_lengths_sum(const Aut& aut);


//! Returns the square root of sum of squares of images lengths differences
LongInteger sqrt_of_sq_dif_of_lengths(const Aut& aut, const Aut& aut1);

//! Returns statistics on the absoulute value of distances of length of images
Statistic<LongInteger> image_length_distance(unsigned int rank, const Aut& e1,
                                              const Aut& e2);



//---------------------------------------------------
// Morphisms enumerators

class AllNielsenGeneratorsEnumerator {
  public:
    AllNielsenGeneratorsEnumerator(unsigned int rank)
      : rank_(rank) {}

    template<typename AutAcceptor>
    void operator()(AutAcceptor& f) {
      Aut::for_each_basic_morphism(rank_, f);
    }

  private:
    unsigned int rank_;
};


class MultiplicationsEnumerator {
  public:
    MultiplicationsEnumerator(unsigned int rank)
      : rank_(rank) {
    }

    template<typename AutAcceptor>
    void operator()(AutAcceptor& f) {
      Aut::for_each_multiplication(rank_, f);
    }

  private:
    unsigned int rank_;
};

class NielsenCompositionGeneratorsEnumerator {
  public:
    NielsenCompositionGeneratorsEnumerator(unsigned int rank)
      : rank_(rank) {}

    template<typename AutAcceptor>
    void operator()(AutAcceptor& f) {
      Aut::for_each_basic_morphism(rank_,
                                   [&] (const Aut& aut) {
        Aut::for_each_basic_morphism(rank_, [&] (const Aut& aut1) {
          f(aut * aut1);
        });
      });
    }

  private:
    unsigned int rank_;
};


//----------------------------------------------------
// Structures holding results

//! Result of minimization
struct MinimizationResult {
  std::vector<Aut> min_morphisms_;
  std::vector<Aut > minimizing_sequence_;
  LongInteger min_value_;

  MinimizationResult() {}

  //! Load from file
  MinimizationResult(std::istream* in);

  //! Saves to stream so it can be recovered later with load
  void save(std::ostream* out) const;
};

//! Call the free reduction for the interval [begin, end) and output to out
template<typename In, typename Out>
void freely_reduce(In begin, In end, Out out) {
  std::transform(begin, end,
                 out,
                 [] (const Aut& aut) {return aut.free_reduction();});
}

template<typename In, typename Func>
Statistic<LongInteger> get_statistic(In begin, In end, Func f) {
  Statistic<LongInteger> stat;
  std::for_each(begin, end,
                [&] (const Aut& aut) {
                  stat(f(aut));
                });
  return stat;
}



extern const char COMMENT_LINE_START;

template<typename T>
void write_comment(std::ostream* out, const T& msg) {
  *out << COMMENT_LINE_START << msg << std::endl;
}

//! Minimization
template<typename Func, typename MinEnum>
MinimizationResult minimize_morphisms(std::ostream* p_log, unsigned int rank,
                                        const std::vector<Aut>& morphisms,
                                        Func f,
                                        MinEnum minimizators_enumerator,
                                        bool logging = false) {
  std::stringstream log_buffer;
  auto comment = [&] (const LongInteger& value, const std::vector<Aut>& morphs) {
    log_buffer << "val=" << value << ",";
    unsigned int max_height = 0;
    unsigned long total_vertices_num = 0;
    std::for_each(morphs.cbegin(), morphs.cend(), [&] (const Aut& aut) {
      auto h = height(aut);
      auto v_num = slp_vertices_num(aut);
      max_height = h > max_height ? h : max_height;
      total_vertices_num += v_num;
    });
    log_buffer << "max_height=" << max_height << ", total vert num=" << total_vertices_num;
    write_comment(p_log, log_buffer.str());
    log_buffer.str(std::string());
  };

  std::size_t num = morphisms.size();
  assert(num > 0);
  MinimizationResult result;
  result.min_morphisms_.reserve(num);
  freely_reduce(morphisms.cbegin(), morphisms.cend(), std::back_inserter(result.min_morphisms_));

  auto min_value = result.min_value_ = f(result.min_morphisms_);
  if (logging) {
    log_buffer << "start ";
    comment(result.min_value_, result.min_morphisms_);
  }

  std::vector<Aut> min_trial;
  Aut minimizing_conjugator;

  auto minimizer = [&] (const Aut& e) {
    std::vector<Aut> trial;
    trial.reserve(num);
    std::transform(result.min_morphisms_.cbegin(), result.min_morphisms_.cend(),
                   std::back_inserter(trial),
                   [&e] (const Aut& aut) {
                      return (e.inverse() * aut * e).free_reduction();
                   });

    auto value = f(trial);
    if (value < min_value) {
      min_value = value;
      min_trial = trial;
      minimizing_conjugator = e;
    }
  };

  while(true) {
    std::vector<Aut> normal_forms;
    normal_forms.reserve(num);
    std::transform(result.min_morphisms_.cbegin(), result.min_morphisms_.cend(),
                   std::back_inserter(normal_forms),
                   [&] (const Aut& aut) {
                      return aut.normal_form();
                   });
    result.min_morphisms_ = normal_forms;
    if (logging) {
      log_buffer << "normal forms ";
      comment(result.min_value_, result.min_morphisms_);
    }

    minimizators_enumerator(minimizer);

    if (result.min_value_ <= min_value) {
      if (logging) {
        write_comment(p_log, "Could not minimize more. Finishing.");
      }
      break;
    }
    if (logging) {
      log_buffer << "minimization step ";
      comment(min_value, min_trial);
    }
    result.min_morphisms_ = min_trial;
    result.min_value_ = min_value;
    result.minimizing_sequence_.push_back(minimizing_conjugator);
  }
  return result;
}


//
struct AutDecomposition {
    Aut aut;
    std::vector<Aut> decomposition;
    std::vector<Aut> conjugators;

    AutDecomposition()
      : aut(), decomposition(), conjugators() {}

    //! Generates a random autmorphism.l
    template<typename RandomAutomorphismGenerator>
    AutDecomposition(unsigned int size, RandomAutomorphismGenerator& random) {
      decomposition.reserve(size);
      for (unsigned int i = 0; i < size; ++i) {
        decomposition.push_back(random());
      }
      aut = Aut::composition(decomposition.begin(), decomposition.end());
    }

    AutDecomposition(std::istream* p_in) {
      //todo
    }

    std::size_t num() const {
      return decomposition.size();
    }

    template<typename It>
    void conjugate_with(It begin, It end) {
      aut.conjugate_with(begin, end);
      std::copy(begin, end,
                std::back_inserter(conjugators));
    }

    AutDecomposition inverse() const {
      AutDecomposition a_dec;
      a_dec.decomposition.reserve(decomposition.size());
      std::for_each(decomposition.rbegin(), decomposition.rend(), [&] (const Aut& a) {
        Aut inverse = a.inverse();
        a_dec.decomposition.push_back(inverse);
        a_dec.aut *= inverse;
      });

      a_dec.aut.conjugate_with(conjugators.begin(), conjugators.end());
      a_dec.conjugators.reserve(conjugators.size());
      std::copy(conjugators.begin(), conjugators.end(),
                std::back_inserter(a_dec.conjugators));

      return a_dec;
    }

    void save(std::ostream* p_out) const {
      //todo
    }
};


//---------------------------------------------------------------------------------------------------------
// experiments


enum GenerateImages {
  generate_images,
  not_generate_images
};

struct Result {
  unsigned int aut_num_;
  unsigned int rank_;
  unsigned int comp_num_;
  unsigned int conj_num_;

  std::vector<AutDecomposition > morphisms_;
  MinimizationResult minimized_morphisms_;
  std::vector<Aut> conjugations_;
  LongInteger conjugation_value_;
  MinimizationResult min_conjugations_;

  std::vector<Aut> conjugation_parts_;
  //conjugator = prod of conjugation parts


  //! Generates experiment data and runs experiment.
  template<typename Generator, typename TargetFunc, typename MinEnumerator>
  Result(unsigned int aut_num, unsigned int rank, unsigned int comp_num, unsigned int conj_num, Generator& generator, TargetFunc target_func, MinEnumerator minimizators_enumerator, std::ostream* p_log)
    : aut_num_(aut_num),
      rank_(rank),
      comp_num_(comp_num),
      conj_num_(conj_num),
      morphisms_(),
      conjugations_(),
      min_conjugations_() {
    conjugation_parts_.reserve(conj_num);
    for (unsigned int i = 0; i < conj_num; ++i) {
      conjugation_parts_.push_back(generator());
    }

    morphisms_.reserve(aut_num);
    std::vector<Aut> v;
    v.reserve(aut_num);
    for (unsigned int i = 0; i < aut_num; ++i) {
      AutDecomposition comp(comp_num, generator);
      morphisms_.push_back(comp);
      v.push_back(comp.aut);
      minimized_morphisms_ = minimize_morphisms(p_log,
                                                rank,
                                               v,
                                               target_func,
                                               minimizators_enumerator, false);
    }


    const std::vector<Aut>& m_morphs = minimized_morphisms_.min_morphisms_;

    std::for_each(m_morphs.begin(), m_morphs.end(), [&] (const Aut& aut) {
      conjugations_.push_back(aut.conjugate_with(conjugation_parts_.begin(), conjugation_parts_.end()));
    });

    conjugation_value_ = target_func(conjugations_);

    min_conjugations_ = minimize_morphisms(p_log,
                                           rank,
                                    conjugations_,
                                    target_func,
                                    minimizators_enumerator, true);

  }

  //! Load from file
  Result(std::istream* in);

  //! Saves to stream so it can be recovered later with load
  void save(std::ostream* out) const;
  //! Prints in more user friendly format but can not be recovered
  /**
   * Uses graphviz representation for graphs.
   */
  void print(std::ostream* out) const;
  /**
   * Prints html representation.
   */
  void print_html(std::ostream* p_html, const std::string& dir, const std::string& filename_prefix, unsigned int exp_index, GenerateImages is_generating_images) const;
};

//! Results of experiment along with target minimization value.
struct Results {
    std::string value_name;
    std::vector<Result> exp_results;
};

//! Reads results from file.
Results read_results(const std::string& filename);

//! Prints morphism with SLP of height less than 3
void print_basic_morphism(std::ostream* p_out, const Aut& e);

//! Prints images of generators under #aut
void print_explicit_images(std::ostream* p_out, const Aut& aut);

void print_values_to_csv_file(const string& filename, const Results& results);


struct indent {
    unsigned int level;

    indent(unsigned int level) : level(level) {}

    friend std::ostream& operator<<(std::ostream& stream, const indent& val) {
      for(int i = 0; i < val.level; ++i) {
          stream << "  ";
      }
      return stream;
  }
};



template<typename In, typename In1>
std::size_t num_of_equal_morphisms(In begin, In end, In1 begin1) {
  std::size_t n(0);
  std::for_each(begin, end, [&] (const Aut& aut) {
    auto aut1 = *begin1;
    ++begin1;
    if (aut == aut1) {
      ++n;
    }
  });
  return n;
}


template<typename Filter>
void print_html(const std::string& dir, const std::string& filenames_prefix, const Results& results, GenerateImages is_generating_images, Filter filter = [] (const Result& r) {return true;}) {
  assert (results.exp_results.size() > 0);
  std::string stylesheet_filename = "stylesheet.css";

  std::ofstream style(dir + stylesheet_filename);

  style << "th {" << std::endl;
  style << indent(1) << "padding: 5px;" << std::endl;
  style << indent(1) << "border-left: 1px solid black;" << std::endl;
  style << indent(1) << "border-bottom: 1px solid black;" << std::endl;
  style << "}" << std::endl;

  style << "th.left {" << std::endl;
  style << indent(1) << "border-left: none;" << std::endl;
  style << "}" << std::endl;

  style << "td {" << std::endl;
  style << indent(1) << "text-align: center;" << std::endl;
  style << "}" << std::endl;


  std::ofstream html(dir + "index.html");
  html << "<!DOCTYPE html>" << std::endl;
  html << "<html>" << std::endl;
  html << "<head>" << std::endl;
  html << "<link type=\"text/css\" rel=\"stylesheet\" href=\"" << stylesheet_filename << "\" />" << std::endl;
  html << "  <title>" << "Results for minimization value " << results.value_name << "</title>" << std::endl;
  html << "</head>" << std::endl;

  html << "<body>" << std::endl;

  const std::size_t tuple_size = results.exp_results.front().aut_num_;

  std::vector<unsigned int> eq_nums(tuple_size + 1, 0);

  //preliminary statistics
  for (const Result& result: results.exp_results) {
    if (filter(result)) {
      const auto& morphs = result.minimized_morphisms_.min_morphisms_;
      const auto& mins = result.min_conjugations_.min_morphisms_;
      auto n = num_of_equal_morphisms(morphs.cbegin(), morphs.cend(), mins.cbegin());
      assert (n >=0 && n <= tuple_size);
      eq_nums[n]++;
    }
  }

  html << "<h3>Number of equal morphisms</h3>" << std::endl;
  html << "<ul>" << std::endl;
  for (std::size_t i = 0; i < eq_nums.size(); ++i) {
    html << indent(1) << "<li>" << i << ": " << eq_nums[i] << "</li>" << std::endl;
  }
  html << "</ul>" << std::endl;


  int n = 0;
  for (const Result& result: results.exp_results) {
    if (filter(result)) {
      result.print_html(&html, dir, filenames_prefix, n++, is_generating_images);
    }
  }

  html << "</body>" << std::endl;

  html << "</html>" << std::endl;

}


void print_all_html(const std::string& dir, const std::string& filenames_prefix, const Results& results, GenerateImages is_generating_images);

void print_not_successful_html(const std::string& dir, const std::string& filenames_prefix, const Results& results, GenerateImages is_generating_images);



// -----------------------------------------------------------------------
// Kicking minimization result out of local minima to get better result.

struct KickResult {
  unsigned int aut_num_;
  unsigned int conj_num_;
  unsigned int rank_;
  std::vector<Aut> conjugations_;
  MinimizationResult min_conjugations_;

  std::vector<Aut> conjugation_parts_;
  //conjugator = prod of conjugation parts


  template<typename Generator, typename TargetFunc, typename MinEnumerator>
  KickResult(const std::vector<Aut>& morphisms, unsigned int rank,  unsigned int conj_num, Generator& generator, TargetFunc target_func, MinEnumerator minimizators_enumerator, std::ostream* p_log)
    : aut_num_(morphisms.size()),
      conj_num_(conj_num),
      rank_(rank)
  {
    conjugation_parts_.reserve(conj_num);
    for (unsigned int i = 0; i < conj_num; ++i) {
      conjugation_parts_.push_back(generator());
    }

    std::for_each(morphisms.begin(), morphisms.end(), [&] (const Aut& aut) {
      conjugations_.push_back(aut.conjugate_with(conjugation_parts_.begin(), conjugation_parts_.end()));
    });

    min_conjugations_ = minimize_morphisms(p_log,
                                           rank,
                                    conjugations_,
                                    target_func,
                                    minimizators_enumerator, true);

  }

  //! Load from file
  KickResult(std::istream* in);

  //! Saves to stream so it can be recovered later with load
  void save(std::ostream* out) const;
  /**
   * Prints html representation.
   */
  void print_html(std::ostream* p_html) const;
};

std::pair<Result, std::vector<KickResult> > read_kick_results(const std::string& filename);

void print_kick_results_to_html(const std::string& dir, const std::string& filenames_prefix, const Result& result, const std::vector<KickResult>& kick_results);

} //namespace fga_csp_attack
} //namespace crag


#endif // CRAG_FREE_GROUPS_FGA_CSP_ATTACK_H
