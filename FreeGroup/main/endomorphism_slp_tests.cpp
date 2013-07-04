
/*
 * endomorphism_slp_tests.cpp
 *
 *  Created on: Mar 19, 2013
 *      Author: pmorar
 */
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>
#include <functional>

#include <cstdlib>

#include <sys/stat.h>
#include <pwd.h>

#include "EndomorphismSLP.h"

typedef std::chrono::high_resolution_clock our_clock;
using namespace crag;


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

typedef EndomorphismSLP<int> Aut;

// Lenght functions for automorphisms.

Statistic<LongInteger> get_endomorphism_images_lengths_stat(const Aut& e) {
  Statistic<LongInteger> length_stat;
  std::map<int, LongInteger> lengths = images_length(e);
  for (const auto& pair: lengths) {
    length_stat(pair.second);
  }
  return length_stat;
}

LongInteger total_images_length(const Aut& em) {
  return get_endomorphism_images_lengths_stat(em).sum();
}

LongInteger max_images_length(const Aut& em) {
  return get_endomorphism_images_lengths_stat(em).max();
}

LongInteger sqrt_of_sq_of_img_lengths_sum(const Aut& aut) {
  LongInteger sum(0);
  std::map<int, LongInteger> lengths = images_length(aut);
  for (const auto& pair: lengths) {
    sum += pair.second * pair.second;
  }
  return sqrt(sum);
}

LongInteger sqrt_of_sq_dif_of_lengths(const Aut& aut, const Aut& aut1) {
  LongInteger sum(0);
  std::map<int, LongInteger> lengths = images_length(aut);
  std::map<int, LongInteger> lengths1 = images_length(aut1);
  std::set<int> parsed_keys;
  for (const auto& pair: lengths) {
    auto key = pair.first;
    auto val = pair.second;
    auto k_v = lengths1.find(key);
    if (k_v != lengths1.end()) {
      auto val1 = k_v->second;
      auto d = val - val1;
      sum += d * d;
    } else {
      sum + val * val;
    }
    parsed_keys.insert(key);
  }

  for (const auto& pair: lengths1) {
    if (parsed_keys.find(pair.first) != parsed_keys.end()) {
      continue;
    }
    auto val = pair.second;
    sum + val * val;
  }

  return sqrt(sum);
}

//! returns statistics on the absoulute value of distances of length of images
Statistic<LongInteger> image_length_distance(unsigned int rank, const Aut& e1,
                                              const Aut& e2) {
  Statistic<LongInteger> length_stat;
  for (int i = 1; i <= rank; ++i) {
    auto img1 = e1.image(i);
    auto img2 = e2.image(i);
    LongInteger d = img1.length() - img2.length();
    length_stat(d > 0 ? d : -d);
  }
  return length_stat;
}



void print_stats(std::ostream* out, const Aut& e) {
  *out << "height=" << height(e) << ", vertices num=" << slp_vertices_num(e)
       << ", image lengths=(" << get_endomorphism_images_lengths_stat(e) << ")";
}

template<typename In>
void print_stats(std::ostream* out, In begin, In end) {
  *out << "{" << std::endl;
  std::for_each(begin, end, [&] (const Aut& aut) {
    print_stats(out, aut);
    *out << std::endl;
  });
  *out << "}" << std::endl;
}


const char COMMENT_LINE_START = '#';
const char* EXPERIMENT_DELIMITER = "#-------------------------";
const char* MORPHISMS_DELIMITER = "***";

template<typename T>
void write_comment(std::ostream* out, const T& msg) {
  *out << COMMENT_LINE_START << msg << std::endl;
}

void skip_comments(std::istream* in) {
  while (in->peek() == COMMENT_LINE_START ||
         in->peek() == '\n') {
    in->ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  }
}



// Minimization


struct MinimizationResult {
  std::vector<Aut> min_morphisms;
  std::vector<Aut > minimizing_sequence;
  LongInteger min_value;

  MinimizationResult() {}

  //! Load from file
  MinimizationResult(std::istream* in);

  //! Saves to stream so it can be recovered later with load
  void save(std::ostream* out) const;
};




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
  result.min_morphisms.reserve(num);
  freely_reduce(morphisms.cbegin(), morphisms.cend(), std::back_inserter(result.min_morphisms));

  auto min_value = result.min_value = f(result.min_morphisms);
  if (logging) {
    log_buffer << "start ";
    comment(result.min_value, result.min_morphisms);
  }

  std::vector<Aut> min_trial;
  Aut minimizing_conjugator;

  auto minimizer = [&] (const Aut& e) {
    std::vector<Aut> trial;
    trial.reserve(num);
    std::transform(result.min_morphisms.cbegin(), result.min_morphisms.cend(),
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
    std::transform(result.min_morphisms.cbegin(), result.min_morphisms.cend(),
                   std::back_inserter(normal_forms),
                   [&] (const Aut& aut) {
                      return aut.normal_form();
                   });
    result.min_morphisms = normal_forms;
    if (logging) {
      log_buffer << "normal forms ";
      comment(result.min_value, result.min_morphisms);
    }

    minimizators_enumerator(minimizer);

    if (result.min_value <= min_value) {
      if (logging) {
        write_comment(p_log, "Could not minimize more. Finishing.");
      }
      break;
    }
    if (logging) {
      log_buffer << "minimization step ";
      comment(min_value, min_trial);
    }
    result.min_morphisms = min_trial;
    result.min_value = min_value;
    result.minimizing_sequence.push_back(minimizing_conjugator);
  }
  return result;
}



template<typename Aut>
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
      a_dec.conjugators.reverse(conjugators.size());
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

void composition_statistics() {
  typedef unsigned int uint;
  std::cout << "Legend: (num iterations, num of composed elements, rank)" << std::endl;
  for (auto rank : {3, 5, 10, 20}) {
    UniformAutomorphismSLPGenerator<int> rnd(rank);
    for (auto size : {100, 1000, 2000}) {
      const uint iterations_num = 100;
      Statistic<unsigned int> height_stat;
      Statistic<unsigned int> vertices_num_stat;

      our_clock::duration time;
      for (int i = 0; i < iterations_num; ++i) {
        auto start_time = our_clock::now();
        auto e = Aut::composition(size, rnd);
        time += our_clock::now() - start_time;
        auto height = crag::height(e);
        auto vertices_num = slp_vertices_num(e);

        height_stat(height);
        vertices_num_stat(vertices_num);
      }
      auto time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(time);
      std::cout << "(iterations=" << iterations_num << ", size=" << size << ", rank=" << rank << "): "
                << time_in_ms.count() << "ms, "
                << time_in_ms.count() / iterations_num << "ms per iteration, "
                << "height " << height_stat << ", "
                << "vertices num " << vertices_num_stat
                << std::endl;
    }
  }
}

void normal_form_statistics(const std::string& filename) {
  std::ofstream out(filename);
  typedef unsigned int uint;
  out << "rank;|e|;vertices_num;height;free_red_vn;free_red_h;jez_nf_vn;jez_ht;jez_of_fr_vn;jez_of_fr_h;fr_of_jez_vn;fr_of_jez_h;" << std::endl;

  auto print_stats = [&out] (const Aut& aut) {
    auto h = height(aut);
    auto vertices_num = slp_vertices_num(aut);
    out << vertices_num << ";" << h << ";";
  };

  for (auto rank : {3, 5}) {
    UniformAutomorphismSLPGenerator<int> rnd(rank);
    for (auto size : {640, 1280}) {
      const uint iterations_num = 10;

      our_clock::duration comp_time;
      our_clock::duration free_red_time;
      our_clock::duration jez_nf_time;
      our_clock::duration jez_of_free_red_time;
      our_clock::duration free_red_of_jez_time;
      for (int i = 0; i < iterations_num; ++i) {
        auto start_time = our_clock::now();
        auto e = Aut::composition(size, rnd);
        comp_time += our_clock::now() - start_time;

        start_time = our_clock::now();
        auto free_red = e.free_reduction();
        free_red_time += our_clock::now() - start_time;

        start_time = our_clock::now();
        auto nf = e.normal_form();
        jez_nf_time += our_clock::now() - start_time;

        start_time = our_clock::now();
        auto nf_of_free_red = Aut();//free_red.normal_form();
        jez_of_free_red_time += our_clock::now() - start_time;

        start_time = our_clock::now();
        auto free_red_of_jez = Aut();//nf.free_reduction();
        free_red_of_jez_time += our_clock::now() - start_time;

        out << rank << ";" << size << ";";
        print_stats(e);
        print_stats(free_red);
        print_stats(nf);
        print_stats(nf_of_free_red);
        print_stats(free_red_of_jez);
        out << std::endl;
      }
      auto comp_time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(comp_time);
      auto free_red_time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(free_red_time);
      auto jez_nf_time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(jez_nf_time);
      auto jez_of_free_red_time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(jez_of_free_red_time);
      auto free_red_of_jez_time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(free_red_of_jez_time);

      std::cout << "(rank=" << rank << ", iterations=" << iterations_num << ", |e|=" << size << std::endl;
      std::cout << ", comp_time=" << comp_time_in_ms.count() << "ms, " <<
                   "free_red_time=" << free_red_time_in_ms.count() << "ms, " <<
                   "jez_nf_time=" << jez_nf_time_in_ms.count() << "ms, " <<
                   "jez_of_free_red_time=" << jez_of_free_red_time_in_ms.count() << "ms, " <<
                   "free_red_of_jez_time=" << free_red_of_jez_time_in_ms.count() << "ms, " << std::endl;
    }
  }
}


enum GenerateImages {
  generate_images,
  not_generate_images
};

struct Result {
  unsigned int aut_num_;
  unsigned int rank_;
  unsigned int comp_num_;
  unsigned int conj_num_;

  std::vector<AutDecomposition<Aut> > morphisms_;
  MinimizationResult minimized_morphisms_;
  std::vector<Aut> conjugations_;
  LongInteger conjugation_value_;
  MinimizationResult min_conjugations_;

  std::vector<Aut> conjugation_parts_;
  //conjugator = prod of conjugation parts


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
      AutDecomposition<Aut> comp(comp_num, generator);
      morphisms_.push_back(comp);
      v.push_back(comp.aut);
      minimized_morphisms_ = minimize_morphisms(p_log,
                                                rank,
                                               v,
                                               target_func,
                                               minimizators_enumerator, false);
    }


    const std::vector<Aut>& m_morphs = minimized_morphisms_.min_morphisms;

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

void print_basic_morphism(std::ostream* p_out, const Aut& e) {
  e.for_each_non_trivial_image([&p_out] (const Aut::symbol_image_pair_type& pair) {
    *p_out << pair.first << " -> ";
    slp::Vertex v = pair.second;
    if (v.height() <= 1)
      *p_out << slp::TerminalVertexTemplate<int>(v).terminal_symbol();
    else {
      auto l = v.left_child();
      auto r = v.right_child();
      *p_out << slp::TerminalVertexTemplate<int>(l).terminal_symbol() << " "
                << slp::TerminalVertexTemplate<int>(r).terminal_symbol();
    }
    *p_out << std::endl;
  });
}



MinimizationResult::MinimizationResult(std::istream* p_in) {
  std::istream& in = *p_in;
  skip_comments(p_in);
  in >> min_value;
  skip_comments(p_in);
  int num;
  in >> num;
  min_morphisms.reserve(num);
  for (unsigned int i = 0; i < num; ++i) {
    skip_comments(p_in);
    min_morphisms.push_back(Aut::load_from(p_in));
  }
  skip_comments(p_in);
  in >> num;

  minimizing_sequence.reserve(num);
  for (unsigned int i = 0; i < num; ++i) {
    skip_comments(p_in);
    minimizing_sequence.push_back(Aut::load_from(p_in));
  }
}

void MinimizationResult::save(std::ostream* p_out) const {
  std::ostream& out = *p_out;
  write_comment(p_out, "min value");
  out << min_value << std::endl;
  write_comment(p_out, "minimized elements");
  write_comment(p_out, "num");
  out << min_morphisms.size() << std::endl;
  write_comment(p_out, "");
  for (const auto& morphism: min_morphisms) {
    morphism.save_to(p_out);
    write_comment(p_out, "");
  }

  write_comment(p_out, "minimizing sequence");
  write_comment(p_out, "num");
  out << minimizing_sequence.size() << std::endl;
  write_comment(p_out, "");
  for (const auto& morphism: minimizing_sequence) {
    morphism.save_to(p_out);
    write_comment(p_out, "");
  }


}

void Result::save(std::ostream* p_out) const {
  std::ostream& out = *p_out;

  out << EXPERIMENT_DELIMITER << std::endl;

  write_comment(p_out, "num of automorphisms");
  out << aut_num_ << std::endl;

  write_comment(p_out, "rank");
  out << rank_ << std::endl;

  write_comment(p_out, "num of composed morphisms");
  out << comp_num_ << std::endl;

  write_comment(p_out, "num of conjugation morphisms");
  out << conj_num_ << std::endl;

  write_comment(p_out, "conjugation value");
  out << conjugation_value_ << std::endl;

  write_comment(p_out, "morphisms");
  for (auto& morph: morphisms_) {
    morph.save(p_out);
    write_comment(p_out, "");
  }

  write_comment(p_out, "minimized morphisms");
  minimized_morphisms_.save(p_out);

  write_comment(p_out, "conjugation morphisms");
  for (const auto& morphism: conjugation_parts_) {
    morphism.save_to(p_out);
    write_comment(p_out, "");
  }

  write_comment(p_out, "conjugations");
  for (const auto& morphism: conjugations_) {
    morphism.save_to(p_out);
    write_comment(p_out, "");
  }

  write_comment(p_out, "minimized conjugations");
  min_conjugations_.save(p_out);
}

void print_explicit_images(std::ostream* p_out, const Aut& aut) {
  std::ostream& out = *p_out;
  aut.for_each_non_trivial_image([&out] (const Aut::symbol_image_pair_type& pair) {
    out << "<div>" << std::endl;
    out << pair.first << " -> ";
    slp::VertexWord<int> image(pair.second);
    std::for_each(image.begin(), image.end(), [&out] (int symbol) {out << symbol;});
    out << std::endl;
    out << "</div>" << std::endl;
  });
}

//template<typename T>
//void Result<T>::print(std::ostream* p_out) const {
//  std::ostream& out = *p_out;

//  out << EXPERIMENT_DELIMITER << std::endl;

//  write_comment(p_out, "rank");
//  out << rank_ << std::endl;

//  write_comment(p_out, "original morphism value");
//  out << min_result_.morphism_value << std::endl;

//  write_comment(p_out, "conjugation value");
//  out << min_result_.conjugation_value << std::endl;

//  write_comment(p_out, "minimized morphism value");
//  out << min_result_.minimized_morphism_value << std::endl;

//  write_comment(p_out, "minimized conjugation value");
//  out << min_result_.minimized_conjugation_value << std::endl;

//  write_comment(p_out, "num of composed morphisms");
//  out << composition_parts_.size() << std::endl;

//  write_comment(p_out, "composed morphisms");
//  for (const auto& morphism: composition_parts_) {
//    print_basic_morphism(p_out, morphism);
//  }

//  write_comment(p_out, "num of conjugation morphisms");
//  out << conjugation_parts_.size() << std::endl;

//  write_comment(p_out, "conjugation morphisms");
//  for (const auto& morphism: conjugation_parts_) {
//    print_basic_morphism(p_out, morphism);
//  }

//  write_comment(p_out, "original morphism");
//  morphisms_.save_graphviz(p_out, "Morphism");

//  write_comment(p_out, "conjugation");
//  conjugations_.save_graphviz(p_out, "Conjugation");

//  write_comment(p_out, "minimized morphism");
//  min_result_.morphism.save_graphviz(p_out, "MinMorphism");

//  write_comment(p_out, "minimized conjugation");
//  min_result_.minimized_conjugation.save_graphviz(p_out, "MinConjugation");

//  write_comment(p_out, "num of morphism minimizing conjugators");
//  out << min_result_.morphism_minimizing_sequence.size() << std::endl;

//  write_comment(p_out, "sequence minimizing morphism");
//  for (const auto& morphism: min_result_.morphism_minimizing_sequence) {
//    print_basic_morphism(p_out, morphism);
//  }

//  write_comment(p_out, "num of conjugation minimizing conjugators");
//  out << min_result_.conjugation_minimizing_sequence.size() << std::endl;

//  write_comment(p_out, "sequence minimizing conjugaion");
//  for (const auto& morphism: min_result_.conjugation_minimizing_sequence) {
//    print_basic_morphism(p_out, morphism);
//  }
//}


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

void Result::print_html(std::ostream* p_html, const std::string& dir, const std::string& aux_filename_prefix, unsigned int exp_index, GenerateImages is_generating_images) const {
  std::ostream& html = *p_html;

  html << indent(1) << "<h3>" << "Experiment " << exp_index << "</h3>" << std::endl;

  html << indent(1) << "<p>" << "Parameters:" << "</p>" << std::endl;
  html << indent(1) << "<ul>" << std::endl;
  html << indent(2) << "<li>tuple size =" << aut_num_ << "</li>" << std::endl;
  html << indent(2) << "<li>rank = " << rank_ << "</li>" << std::endl;
  html << indent(2) << "<li>|base| = " << comp_num_ << "</li>" << std::endl;
  html << indent(2) << "<li>|conjugator| =" << conj_num_ << "</li>" << std::endl;
  html << indent(2) << "<li>|minimizing conjugator| = " << min_conjugations_.minimizing_sequence.size() << "</li>" << std::endl;
  html << indent(1) << "</ul>" << std::endl;

  html << indent(1) << "<p>" << "Values:" << "</p>" << std::endl;
  html << indent(1) << "<ul>" << std::endl;
  html << indent(2) << "<li>morphism (minimized) = " << minimized_morphisms_.min_value << "</li>" << std::endl;
  html << indent(2) << "<li>conjugation = " << conjugation_value_ << "</li>" << std::endl;
  html << indent(2) << "<li>minimized conjugation = " << min_conjugations_.min_value << "</li>" << std::endl;
  html << indent(1) << "</ul>" << std::endl;


  html << indent(1) << "<table>" << std::endl;

  html << indent(2) << "<thead>" << std::endl;

//  html << indent(3) << "<tr>" << std::endl;
//  html << indent(4) << "<th class=\"left\" colspan=\"2\">Sample genrating data</th>" << std::endl;
//  html << indent(4) << "<th colspan=\"2\">Minimizing conjugators</th>" << std::endl;
//  html << indent(3) << "</tr>" << std::endl;

  html << indent(3) << "<tr>" << std::endl;
  html << indent(4) << "<th>Conjugators</th>" << std::endl;
  html << indent(4) << "<th>Minimizing conjugators</th>" << std::endl;
  html << indent(3) << "</tr>" << std::endl;

  html << indent(2) << "</thead>" << std::endl;

  html << indent(2) << "<tbody>" << std::endl;
  auto conj_iterator = conjugation_parts_.cbegin();
  auto min_conj_iterator = min_conjugations_.minimizing_sequence.cbegin();

  while (true) {
    bool all_done = true;

    if(conj_iterator != conjugation_parts_.cend()) {
      all_done = false;
      html << indent(4) << "<td>";
      print_basic_morphism(p_html, *conj_iterator);
      html << "</td>" << std::endl;
      ++conj_iterator;
    } else {
      html << "<td></td>" << std::endl;
    }

    if(min_conj_iterator != min_conjugations_.minimizing_sequence.cend()) {
      all_done = false;
      html << indent(4) << "<td>";
      print_basic_morphism(p_html, *min_conj_iterator);
      html << "</td>" << std::endl;
      ++min_conj_iterator;
    } else {
      html << "<td></td>" << std::endl;
    }

    html << indent(3) << "</tr>" << std::endl;
    if (all_done)
      break;
  }

  html << indent(2) << "</tbody>" << std::endl;

  html << indent(1) << "</table>" << std::endl;


  html << indent(1) << "<p>Conjugation morphism: " << std::endl;

  Aut conjugation = Aut::composition(conjugation_parts_.cbegin(), conjugation_parts_.cend());
  print_explicit_images(p_html, conjugation);

  html << indent(1) << "</p>" << std::endl;

  html << indent(1) << "<p>Minimization morphism inverse: " << std::endl;

  Aut minimizator = Aut::composition(min_conjugations_.minimizing_sequence.cbegin(), min_conjugations_.minimizing_sequence.cend());
  print_explicit_images(p_html, minimizator);

  html << indent(1) << "</p>" << std::endl;

//  if (is_generating_images == generate_images) {
//    auto generate_morphism_description = [&] (const Aut& e, const std::string& name, const std::string& description) {
//      std::stringstream s;
//      s << aux_filename_prefix << "_" << name << "_" << exp_index;
//      std::string filename = s.str();
//      std::string description_filename(filename + ".gv");
//      std::string image_filename(filename + ".gif");

//      html << indent(1) << "<p>" << description << "</p>";
//      html << indent(1) << "<img src=\"" << image_filename << "\"/>" << std::endl;

//      std::ofstream description_file(dir + description_filename);
//      e.save_graphviz(&description_file, name);
//      s.str("");
//      s << "dot -Tgif " << dir << description_filename << " -o " << dir << image_filename;
//      std::system(s.str().c_str());
//    };

////    generate_morphism_description(morphisms_, "morphism", "Original morphism.");
////    generate_morphism_description(conjugations_, "conjugation", "Conjugation.");
////    generate_morphism_description(min_result_.morphism, "min_morphism", "Minimized morphism");
////    generate_morphism_description(min_result_.minimized_conjugation, "min_conjugation", "Minimized conjugation.");
//  }
}

std::istream& operator>>(std::istream& in, LongInteger& n) {
  std::string s;
  std::getline(in, s);
  n = LongInteger(s);
}

Result::Result(std::istream* p_in) {
  std::istream& in = *p_in;

  skip_comments(p_in);
  in >> aut_num_;
  skip_comments(p_in);
  in >> rank_;
  skip_comments(p_in);
  in >> comp_num_;
  skip_comments(p_in);
  in >> conj_num_;
  skip_comments(p_in);
  in >> conjugation_value_;

  morphisms_.reserve(aut_num_);
  for (unsigned int i = 0; i < aut_num_; ++i) {
    morphisms_.push_back(AutDecomposition<Aut>(p_in));
  }

  minimized_morphisms_ = MinimizationResult(p_in);

  conjugation_parts_.reserve(conj_num_);
  for (unsigned int i = 0; i < conj_num_; ++i) {
    skip_comments(p_in);
    conjugation_parts_.push_back(Aut::load_from(p_in));
  }

  conjugations_.reserve(aut_num_);
  for (unsigned int i = 0; i < aut_num_; ++i) {
    skip_comments(p_in);
    conjugations_.push_back(Aut::load_from(p_in));
  }
  min_conjugations_ = MinimizationResult(p_in);
}


//std::ostream& operator<<(std::ostream& out, const MinimizationResult& values) {
//  return out << values.min_value << ";" << values.conjugation_value << ";"
//                << values.minimized_morphism_value << ";" << values.minimized_conjugation_value << ";";
//}

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



void conjugation_length_based_attack_statistics(const std::string& filename, unsigned int aut_num, unsigned int rank, unsigned int comp_num, unsigned int conj_num, unsigned int iterations_num) {
  std::ofstream out(filename);
  write_comment(&out, "Legnth-base attack to Conjugation Search Problem for Automorphisms of Free Group");
  write_comment(&out, "minimization value name");
  out << "sum of total images length" << std::endl;

  auto target_func = [] (const std::vector<Aut>& morphs) {
    return get_statistic(morphs.cbegin(), morphs.cend(), total_images_length).sum();
  };


  std::cout << "start with rank=" << rank
            << ", comp_num=" << comp_num
            << ", conj_num=" << conj_num
            << ", it_num=" << iterations_num
            << std::endl;


  AllNielsenGeneratorsEnumerator enumerator(rank);
  UniformAutomorphismSLPGenerator<int> rnd(rank);

  int eq_num = 0;
  auto start_time = our_clock::now();
  for (int i = 0; i < iterations_num; ++i) {
    std::cout << "iteration " << i << std::endl;

    auto iteration_start_time = our_clock::now();

    Result result(aut_num, rank, comp_num, conj_num, rnd, target_func, enumerator, &out);

    result.save(&out);

    auto iteration_time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(our_clock::now() - iteration_start_time);
    write_comment(&out, "time");
    out << COMMENT_LINE_START << iteration_time_in_ms.count() <<  "ms" << std::endl;
  }
  auto time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(our_clock::now() - start_time);
  std::stringstream s;
  s << "(iterations num=" << iterations_num << ",comp num =" << comp_num << ",conj num =" << conj_num << ",rank=" << rank << "): "
    << time_in_ms.count() << "ms, "
    << time_in_ms.count() / iterations_num << "ms per iteration, "
    << std::endl;
  std::string str = s.str();
  std::cout << str;
  write_comment(&out, str);

  std::cout << "eq num " << eq_num << " out of " << iterations_num << std::endl;
}


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

void KickResult::save(std::ostream* p_out) const {
  std::ostream& out = *p_out;

  out << EXPERIMENT_DELIMITER << std::endl;

  write_comment(p_out, "num of automorphisms");
  out << aut_num_ << std::endl;

  write_comment(p_out, "num of conjugation morphisms");
  out << conj_num_ << std::endl;

  write_comment(p_out, "rank");
  out << rank_ << std::endl;

  write_comment(p_out, "conjugation morphisms");
  for (const auto& morphism: conjugation_parts_) {
    morphism.save_to(p_out);
    write_comment(p_out, "");
  }

  write_comment(p_out, "conjugations");
  for (const auto& morphism: conjugations_) {
    morphism.save_to(p_out);
    write_comment(p_out, "");
  }

  write_comment(p_out, "minimized conjugations");
  min_conjugations_.save(p_out);
}


void KickResult::print_html(std::ostream* p_html) const {
  std::ostream& html = *p_html;

  html << indent(1) << "<h3>" << "Kick attempt</h3>" << std::endl;

  html << indent(1) << "<p>" << "Parameters:" << "</p>" << std::endl;
  html << indent(1) << "<ul>" << std::endl;
  html << indent(2) << "<li>tuple size =" << aut_num_ << "</li>" << std::endl;
  html << indent(2) << "<li>|conjugator| =" << conj_num_ << "</li>" << std::endl;
  html << indent(2) << "<li>|minimizing conjugator| = " << min_conjugations_.minimizing_sequence.size() << "</li>" << std::endl;
  html << indent(1) << "</ul>" << std::endl;

  html << indent(1) << "<p>minimized value = " << min_conjugations_.min_value << "</p>" << std::endl;

  html << indent(1) << "<table>" << std::endl;

  html << indent(2) << "<thead>" << std::endl;

//  html << indent(3) << "<tr>" << std::endl;
//  html << indent(4) << "<th class=\"left\" colspan=\"2\">Sample genrating data</th>" << std::endl;
//  html << indent(4) << "<th colspan=\"2\">Minimizing conjugators</th>" << std::endl;
//  html << indent(3) << "</tr>" << std::endl;

  html << indent(3) << "<tr>" << std::endl;
  html << indent(4) << "<th>Kick conjugators</th>" << std::endl;
  html << indent(4) << "<th>Minimizing conjugators</th>" << std::endl;
  html << indent(3) << "</tr>" << std::endl;

  html << indent(2) << "</thead>" << std::endl;

  html << indent(2) << "<tbody>" << std::endl;
  auto conj_iterator = conjugation_parts_.cbegin();
  auto min_conj_iterator = min_conjugations_.minimizing_sequence.cbegin();

  while (true) {
    bool all_done = true;

    if(conj_iterator != conjugation_parts_.cend()) {
      all_done = false;
      html << indent(4) << "<td>";
      print_basic_morphism(p_html, *conj_iterator);
      html << "</td>" << std::endl;
      ++conj_iterator;
    } else {
      html << "<td></td>" << std::endl;
    }

    if(min_conj_iterator != min_conjugations_.minimizing_sequence.cend()) {
      all_done = false;
      html << indent(4) << "<td>";
      print_basic_morphism(p_html, *min_conj_iterator);
      html << "</td>" << std::endl;
      ++min_conj_iterator;
    } else {
      html << "<td></td>" << std::endl;
    }

    html << indent(3) << "</tr>" << std::endl;
    if (all_done)
      break;
  }

  html << indent(2) << "</tbody>" << std::endl;

  html << indent(1) << "</table>" << std::endl;


//  html << indent(1) << "<p>Conjugation morphism: " << std::endl;

//  Aut conjugation = Aut::composition(conjugation_parts_.cbegin(), conjugation_parts_.cend());
//  print_explicit_images(p_html, conjugation);

//  html << indent(1) << "</p>" << std::endl;

//  html << indent(1) << "<p>Minimization morphism inverse: " << std::endl;

//  Aut minimizator = Aut::composition(min_conjugations_.minimizing_sequence.cbegin(), min_conjugations_.minimizing_sequence.cend());
//  print_explicit_images(p_html, minimizator);

//  html << indent(1) << "</p>" << std::endl;

}


KickResult::KickResult(std::istream* p_in) {
  std::istream& in = *p_in;

  skip_comments(p_in);
  in >> aut_num_;
  skip_comments(p_in);
  in >> conj_num_;
  skip_comments(p_in);
  in >> rank_;

  conjugation_parts_.reserve(conj_num_);
  for (unsigned int i = 0; i < conj_num_; ++i) {
    skip_comments(p_in);
    conjugation_parts_.push_back(Aut::load_from(p_in));
  }

  conjugations_.reserve(aut_num_);
  for (unsigned int i = 0; i < aut_num_; ++i) {
    skip_comments(p_in);
    conjugations_.push_back(Aut::load_from(p_in));
  }
  min_conjugations_ = MinimizationResult(p_in);
}



void conjugation_lba_kick_attempt(const std::string& filename, const Result& result, unsigned int conj_num, unsigned int iterations_num) {
  std::ofstream out(filename);
  write_comment(&out, "Attempt to kick out of local minima");
  write_comment(&out, "original result");
  result.save(&out);

  auto target_func = [] (const std::vector<Aut>& morphs) {
    return get_statistic(morphs.cbegin(), morphs.cend(), total_images_length).sum();
  };

  auto rank = result.rank_;
  auto min_conjs_ = result.min_conjugations_.min_morphisms;

  AllNielsenGeneratorsEnumerator enumerator(rank);
  UniformAutomorphismSLPGenerator<int> rnd(rank);

  auto start_time = our_clock::now();
  for (int i = 0; i < iterations_num; ++i) {
    std::cout << "iteration " << i << std::endl;

    auto iteration_start_time = our_clock::now();



    KickResult k_res(min_conjs_, rank, conj_num, rnd, target_func, enumerator, &out);

    k_res.save(&out);

    auto iteration_time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(our_clock::now() - iteration_start_time);
    write_comment(&out, "time");
    out << COMMENT_LINE_START << iteration_time_in_ms.count() <<  "ms" << std::endl;
  }
  auto time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(our_clock::now() - start_time);
  std::stringstream s;
  s << "(iterations num=" << iterations_num << ",conj num =" << conj_num << ",rank=" << rank << "): "
    << time_in_ms.count() << "ms, "
    << time_in_ms.count() / iterations_num << "ms per iteration, "
    << std::endl;
  std::string str = s.str();
  std::cout << str;
  write_comment(&out, str);
}


//void lba_success_precentage(const std::string& filename) {
//  std::ofstream out(filename);
//  write_comment(&out, "Legnth-base attack to Conjugation Search Problem for Automorphisms of Free Group");
//  write_comment(&out, "minimization value name");
//  out << "total image legnth" << std::endl;


//  auto target_func = [] (const std::vector<Aut>& p_morphs) -> LongInteger {
//    return 1;
//  };

//  typedef unsigned int uint;
//  for (auto rank : {3, 5}) {
//    std::cout << "rank=" << rank << std::endl;
//    AllNielsenGeneratorsEnumerator enumerator(rank);
//    UniformAutomorphismSLPGenerator<int> rnd(rank);
//    for (uint size: {1, 2, 3, 4, 5}) {
//      std::cout << "num of composed automorphisms=" << size << std::endl;
//      for (auto conj_num: {size}) {
//        std::cout << "num of conjugators=" << conj_num << std::endl;
//        const uint iterations_num = 5;

//        auto start_time = our_clock::now();
//        unsigned int success_num = 0;
//        for (unsigned int i = 0; i < iterations_num; ++i) {
////          if (i % 10 == 0)
////            std::cout << "Iteration " << i << std::endl;

//          auto iteration_start_time = our_clock::now();

//          Result<LongInteger> result(rank, size, conj_num, rnd, target_func, enumerator);


////          if (result.min_result_.minimized_morphism_value == result.min_result_.minimized_conjugation_value &&
////              result.min_result_.morphism == result.min_result_.minimized_conjugation)
////            ++success_num; TODO should fix this
////          result.save(&out);
//          auto iteration_time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(our_clock::now() - iteration_start_time);
//          write_comment(&out, "time");
//          out << COMMENT_LINE_START << iteration_time_in_ms.count() <<  "ms" << std::endl;
//        }
//        auto time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(our_clock::now() - start_time);
//        std::cout << "(iterations num=" << iterations_num << ",size=" << size << ",conjug num=" << conj_num << ",rank=" << rank << "): "
//                  << time_in_ms.count() << "ms, "
//                  << time_in_ms.count() / iterations_num << "ms per iteration, "
//                  << "success num=" << success_num
//                  << std::endl;
//      }
//    }
//  }
//}

void skip_line(std::istream* in) {
  in->ignore(std::numeric_limits<std::streamsize>::max(),'\n');
}

struct Results {
    std::string value_name;
    std::vector<Result> exp_results;
};


void print_values_to_csv_file(const string& filename, const Results& results) {
  std::ofstream csv_out(filename);
  csv_out << "val;conj_val;min_val" << std::endl;

  std::for_each(results.exp_results.cbegin(), results.exp_results.cend(), [&csv_out] (const Result& result) {
    csv_out << result.minimized_morphisms_.min_value << ";" <<
      result.conjugation_value_ << ";" <<
      result.min_conjugations_.min_value << std::endl;
  });
}

Results read_results(const std::string& filename) {
  std::ifstream in(filename);
  in.exceptions (std::ios::eofbit | std::ios::failbit |
                 std::ios::badbit);
  skip_comments(&in);

  Results results;
  std::getline(in, results.value_name);


  try {
    while (in && !in.eof()) {
      results.exp_results.push_back(Result(&in));
    }
  } catch(std::ifstream::failure e) {
    //when we can not read anymore we fall off here.
  }

  std::cout << results.exp_results.size() << " of experiments read." << std::endl;

  return results;
}


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
      const auto& morphs = result.minimized_morphisms_.min_morphisms;
      const auto& mins = result.min_conjugations_.min_morphisms;
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


void print_all_html(const std::string& dir, const std::string& filenames_prefix, const Results& results, GenerateImages is_generating_images) {
  print_html(dir, filenames_prefix, results, is_generating_images, [] (const Result& r) {return true;});
}

void print_not_successful_html(const std::string& dir, const std::string& filenames_prefix, const Results& results, GenerateImages is_generating_images) {
  print_html(dir, filenames_prefix, results, is_generating_images, [] (const Result& r) {
    return r.minimized_morphisms_.min_value != r.min_conjugations_.min_value;//TODO equalities
  });
}


void process_length_base_attack() {
  //parameters
  typedef unsigned int uint;
  const uint aut_num = 1;
  const uint rank = 3;
  const uint comp_num = 5;
  const uint conj_num = 5;
  const uint iterations_num = 5;

  //filenames and dirs
  auto myuid = getuid();
  auto mypasswd = getpwuid(myuid);
  std::string dir(mypasswd->pw_dir);
  dir += "/Documents/exp_res1/";

  std::stringstream s;
  s << "lba_total_sum_a" << aut_num << "comp" << comp_num << "conj" << conj_num << "it" << iterations_num;
  std::string name(s.str());
  std::string filename = dir + name + ".txt";
  std::string csv_filename = dir + name + ".csv";
  std::string html_dir = dir + name + "/";

  //work
  conjugation_length_based_attack_statistics(filename, aut_num, rank, comp_num, conj_num, iterations_num);
  auto results = read_results(filename);

  mode_t process_mask = umask(0);
  int result_code = mkdir(html_dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  if (!result_code && errno == EEXIST) {
    std::cerr << "can not create dir!" << std::endl;
  }
  umask(process_mask);

  print_all_html(html_dir, "lba", results, not_generate_images);
  print_values_to_csv_file(csv_filename, results);
}

std::pair<Result, std::vector<KickResult> > read_kick_results(const std::string& filename) {
  std::ifstream in(filename);
  in.exceptions (std::ios::eofbit | std::ios::failbit |
                 std::ios::badbit);
  skip_comments(&in);

  auto pair = std::make_pair(Result(&in), std::vector<KickResult>());

  try {
    while (in && !in.eof()) {
      pair.second.push_back(KickResult(&in));
    }
  } catch(std::ifstream::failure e) {
    //when we can not read anymore we fall off here.
  }

  std::cout << pair.second.size() << " of experiments read." << std::endl;

  return pair;
}

void print_kick_results_to_html(const std::string& dir, const std::string& filenames_prefix, const Result& result, const std::vector<KickResult>& kick_results) {
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
  html << "  <title>" << "Kick results " << "</title>" << std::endl;
  html << "</head>" << std::endl;

  html << "<body>" << std::endl;

  const std::size_t tuple_size = result.aut_num_;

  std::vector<unsigned int> eq_nums(tuple_size + 1, 0);

//  //preliminary statistics
//  for (const Result& result: results.exp_results) {
//    if (filter(result)) {
//      const auto& morphs = result.minimized_morphisms_.min_morphisms;
//      const auto& mins = result.min_conjugations_.min_morphisms;
//      auto n = num_of_equal_morphisms(morphs.cbegin(), morphs.cend(), mins.cbegin());
//      assert (n >=0 && n <= tuple_size);
//      eq_nums[n]++;
//    }
//  }

//  html << "<h3>Number of equal morphisms</h3>" << std::endl;
//  html << "<ul>" << std::endl;
//  for (std::size_t i = 0; i < eq_nums.size(); ++i) {
//    html << indent(1) << "<li>" << i << ": " << eq_nums[i] << "</li>" << std::endl;
//  }
//  html << "</ul>" << std::endl;

  result.print_html(&html, dir, filenames_prefix, 0, GenerateImages::not_generate_images);


  int n = 0;
  for (const KickResult& res: kick_results) {
    res.print_html(&html);
  }

  html << "</body>" << std::endl;

  html << "</html>" << std::endl;

}

void kick_attempt() {
  const uint conj_num = 10;
  const uint iterations_num = 10;

  std::cout << "Starting kick attempt..." << std::endl;

  std::ifstream in("kick_sample.txt");
  Result res(&in);
  std::cout << "Initial sample read." << std::endl;


  //filenames and dirs
  auto myuid = getuid();
  auto mypasswd = getpwuid(myuid);
  std::string dir(mypasswd->pw_dir);
  dir += "/Documents/exp_res1/";

  std::stringstream s;
  s << "kick_attempt_total_sum_conj" << conj_num << "it" << iterations_num;
  std::string name(s.str());
  std::string filename = dir + name + ".txt";
  std::string html_dir = dir + name + "/";

  //work
  conjugation_lba_kick_attempt(filename, res, conj_num, iterations_num);

  std::cout << "Finished attack." << std::endl;

  auto results_pair = read_kick_results(filename);

  std::cout << "Log file read." << std::endl;

  mode_t process_mask = umask(0);
  int result_code = mkdir(html_dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
  if (!result_code && errno == EEXIST) {
    std::cerr << "can not create dir!" << std::endl;
  }
  umask(process_mask);

  print_kick_results_to_html(html_dir, "lba_kick", results_pair.first, results_pair.second);
}

int main() {
  kick_attempt();
//  process_length_base_attack();
//  auto myuid = getuid();
//  auto mypasswd = getpwuid(myuid);
//  std::string dir(mypasswd->pw_dir);
//  dir += "/Documents/exp_results/";

//  normal_form_statistics(dir + "normal_form_stat_large.csv");
}
