
/*
 * endomorphism_slp_tests.cpp
 *
 *  Created on: Mar 19, 2013
 *      Author: pmorar
 */
#include <iostream>
#include <fstream>
#include <chrono>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/mean.hpp>
#include <boost/accumulators/statistics/min.hpp>
#include <boost/accumulators/statistics/max.hpp>
#include <boost/accumulators/statistics/moment.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>

#include <cstdlib>

#include "EndomorphismSLP.h"

typedef std::chrono::high_resolution_clock our_clock;
using namespace crag;
using namespace boost::accumulators;


//----------------------------------------------------------------------------------------
// Statistic accumulators

//! Our statistics accumulator for non LongInteger. If needed add extra tags for more statistics.
template<typename T>
struct Statistic: public accumulator_set<T, stats<tag::min, tag::max, tag::mean>> {
};

//! Specialization for LongInteger. We need this because accumulator does not work correctly with mean and LongInteger
template<>
struct Statistic<LongInteger>: public accumulator_set<LongInteger, stats<tag::count, tag::sum, tag::max>> {
    void operator()(const LongInteger& val) {
      accumulator_set<LongInteger, stats<tag::count, tag::sum, tag::max>>::operator()(val);
      if (count(*this) == 1 || val < min_) {
        min_ = val;
      }
    }

    LongInteger min_;
};

template<typename T>
T minimum(const Statistic<T>& stat) {
  return min(stat);
}

template<>
LongInteger minimum(const Statistic<LongInteger>& stat) {
 return stat.min_;
}

template<typename T>
LongInteger average(const Statistic<T>& stat) {
  return mean(stat);
}

template<>
LongInteger average(const Statistic<LongInteger>& stat) {
  return sum(stat) / count(stat);
}

template <typename T>
std::ostream& operator<<(std::ostream& out, const Statistic<T>& stat) {
  return out << "(avg = " << average(stat) << ", "
                 << "min = " << minimum(stat) << ", "
                 << "max = " << max(stat) << ")";
}

//----------------------------------------------------
// Helper methods

Statistic<LongInteger> get_endomorphism_images_lengths_stat(const EndomorphismSLP<int>& e) {
  Statistic<LongInteger> length_stat;
  auto length_calculator = [&] (const EndomorphismSLP<int>::symbol_image_pair_type& pair) {
    slp::Vertex v = pair.second;
    length_stat(v.length());
  };
  e.for_each_non_trivial_image(length_calculator);
  return length_stat;
}

LongInteger total_images_length(const EndomorphismSLP<int>& em) {
  return sum(get_endomorphism_images_lengths_stat(em));
}

LongInteger _images_length(const EndomorphismSLP<int>& em) {
  return max(get_endomorphism_images_lengths_stat(em));
}

//! returns statistics on the absoulute value of distances of length of images
Statistic<LongInteger> image_length_distance(unsigned int rank, const EndomorphismSLP<int>& e1,
                                              const EndomorphismSLP<int>& e2) {
  Statistic<LongInteger> length_stat;
  for (int i = 1; i <= rank; ++i) {
    auto img1 = e1.image(i);
    auto img2 = e2.image(i);
    LongInteger d = img1.length() - img2.length();
    length_stat(d > 0 ? d : -d);
  }
  return length_stat;
}

void print_stats(std::ostream* out, const EndomorphismSLP<int>& e) {
  *out << "height=" << height(e) << ", vertices num=" << slp_vertices_num(e)
       << ", image lengths=(" << get_endomorphism_images_lengths_stat(e) << ")";
}




template<typename Func>
std::pair<EndomorphismSLP<int>, std::vector<EndomorphismSLP<int>>> minimize_morphism(unsigned int rank, const EndomorphismSLP<int>& e, Func f, bool logging = false) {
  auto morphism = e.free_reduction();
  auto start_value = f(morphism);
  if (logging) {
    std::cout << "start value = " << start_value << ", ";
    print_stats(&std::cout, morphism);
    std::cout << std::endl;
  }
  std::vector<EndomorphismSLP<int>> minimization_sequence;
  while(true) {
    auto min_value = start_value;
    EndomorphismSLP<int> min_trial = EndomorphismSLP<int>::identity();
    EndomorphismSLP<int> minimizing_conjugator = EndomorphismSLP<int>::identity();

    EndomorphismSLP<int>::for_each_basic_morphism(rank,
                                       [&] (const EndomorphismSLP<int>& e) {
      auto trial = e * morphism * e.inverse();
      auto reduced_trial = trial.free_reduction();
      auto value = f(reduced_trial);
      if (value < min_value) {
        min_value = value;
        min_trial = reduced_trial;
        minimizing_conjugator = e;
      }
    });

    if (start_value <= min_value) {
      if (logging)
        std::cout << "Could not decrease target value. Finishing procedure." << std::endl;
      break;
    }
    if (logging) {
      std::cout << "Success: value=" << min_value << ", ";
      print_stats(&std::cout, min_trial);
      std::cout << std::endl;
    }
    morphism = min_trial;
    start_value = min_value;
    minimization_sequence.push_back(minimizing_conjugator);
  }
  return std::make_pair(morphism, minimization_sequence);
}

template<typename T>
struct TargetValues {
  T morphism_value;
  T conjugation_value;
  T minimized_morphism_value;
  T minimized_conjugation_value;
  EndomorphismSLP<int> minimized_morphism;
  EndomorphismSLP<int> minimized_conjugation;

  std::vector<EndomorphismSLP<int> > morphism_minimizing_sequence;
  std::vector<EndomorphismSLP<int> > conjugation_minimizing_sequence;

  TargetValues() : morphism_value(), conjugation_value(),
    minimized_morphism_value(), minimized_conjugation_value(),
    minimized_morphism(EndomorphismSLP<int>::identity()),
    minimized_conjugation(EndomorphismSLP<int>::identity()),
    morphism_minimizing_sequence(),
    conjugation_minimizing_sequence() {}

  TargetValues(T morphism_val, T conjug_val, T min_morphism_val, T min_conjug_val,
               const EndomorphismSLP<int>& minimized_morph, const EndomorphismSLP<int>& minimized_conj,
               const std::vector<EndomorphismSLP<int> >& morph_minimizing_sequence,
               const std::vector<EndomorphismSLP<int> >& conj_minimizing_sequence)
    : morphism_value(morphism_val), conjugation_value(conjug_val),
      minimized_morphism_value(min_morphism_val), minimized_conjugation_value(min_conjug_val),
      minimized_morphism(minimized_morph), minimized_conjugation(minimized_conj),
      morphism_minimizing_sequence(morph_minimizing_sequence),
      conjugation_minimizing_sequence(conj_minimizing_sequence) {}

  TargetValues(T&& morphism_val, T&& conjug_val, T&& min_morphism_val, T&& min_conjug_val,
               EndomorphismSLP<int>&& minimized_morph, EndomorphismSLP<int>&& minimized_conj,
               std::vector<EndomorphismSLP<int> >&& morph_minimizing_sequence,
               std::vector<EndomorphismSLP<int> >&& conj_minimizing_sequence)
    : morphism_value(morphism_val), conjugation_value(conjug_val),
      minimized_morphism_value(min_morphism_val), minimized_conjugation_value(min_conjug_val),
      minimized_morphism(minimized_morph), minimized_conjugation(minimized_conj),
      morphism_minimizing_sequence(morph_minimizing_sequence),
      conjugation_minimizing_sequence(conj_minimizing_sequence) {}
};

template<typename TargetFunc, typename TargetValueType = LongInteger>
TargetValues<TargetValueType> conjugation_function_minimization_attack(unsigned int rank, const EndomorphismSLP<int>& e,
                                     const EndomorphismSLP<int>& e_conjugated,
                                     TargetFunc target_func,
                                     bool detailed_logging = true) {
  auto reduced_e = e.free_reduction();
  auto reduced_conj = e_conjugated.free_reduction();
  if (detailed_logging)
    std::cout << "Minimizing conjugation..." << std::endl;

  auto min_conj_result = minimize_morphism(rank, reduced_conj, target_func, detailed_logging);
  auto min_conj = min_conj_result.first;
  auto min_conj_seq = min_conj_result.second;

  if (detailed_logging)
    std::cout << "Minimizing e itself..." << std::endl;
  auto min_morph_result = minimize_morphism(rank, reduced_e, target_func, detailed_logging);
  auto min_e = min_morph_result.first;
  auto min_e_seq = min_morph_result.second;

  if (detailed_logging)
    std::cout << "Comparing the conjugation, e, and minimized versions." << std::endl;

  auto reduced_min_conj = min_conj.free_reduction();
  auto reduced_min_e = min_e.free_reduction();

  auto e_val = target_func(reduced_e);
  auto conj_val = target_func(reduced_conj);
  auto min_e_val = target_func(min_e);
  auto min_val = target_func(min_conj);


  if (detailed_logging) {
    std::cout << "e :        ";
    print_stats(&std::cout, e);
    std::cout << ", val=" << e_val << std::endl;

    std::cout << "conj :     ";
    print_stats(&std::cout, reduced_conj);
    std::cout << ", val=" << conj_val << std::endl;

    std::cout << "min e:     ";
    print_stats(&std::cout, reduced_min_e);
    std::cout << ", val=" << min_e_val << std::endl;

    std::cout << "min conj : ";
    print_stats(&std::cout, reduced_min_conj);
    std::cout << ", val=" << min_val << std::endl;
  }



//  slp::MatchingTable mt;
  if (detailed_logging) {
    std::cout << "Images of reduced morphisms" << std::endl;
    for (int i = 1; i <= rank; ++i) {
      std::cout << "image of " << i << std::endl;
      auto img_e = reduced_min_e.image_word(i);
      auto img_min = reduced_min_conj.image_word(i);
      std::cout << "e min :" << img_e << std::endl;
      std::cout << "min   :" << img_min << std::endl;
    }
  }

  return TargetValues<TargetValueType>(e_val, conj_val, min_e_val, min_val, reduced_min_e, reduced_min_conj, min_e_seq, min_conj_seq);
}

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
        auto e = EndomorphismSLP<int>::composition(size, rnd);
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

template<typename T>
struct Result {
  unsigned int rank_;
  EndomorphismSLP<int> morphism_;
  EndomorphismSLP<int> conjugation_;
  TargetValues<T> values_;

  std::vector<EndomorphismSLP<int> > composition_parts_;
  std::vector<EndomorphismSLP<int> > conjugation_parts_;

  Result()
    : rank_(),
      morphism_(EndomorphismSLP<int>::identity()),
      conjugation_(EndomorphismSLP<int>::identity()),
      values_(),
      composition_parts_(),
      conjugation_parts_() {}

  template<typename F>
  Result(unsigned int rank, unsigned int comp_num, unsigned int conj_num, F& f)
    : rank_(rank),
      morphism_(EndomorphismSLP<int>::identity()),
      conjugation_(EndomorphismSLP<int>::identity()),
      values_() {
    composition_parts_.reserve(comp_num);
    conjugation_parts_.reserve(conj_num);
    for (unsigned int i = 0; i < comp_num; ++i)
      composition_parts_.push_back(f());
    for (unsigned int i = 0; i < conj_num; ++i)
      conjugation_parts_.push_back(f());
    morphism_.compose_with(composition_parts_.begin(), composition_parts_.end());
    conjugation_ = morphism_.conjugate_with(conjugation_parts_.begin(), conjugation_parts_.end());
  }

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
  void print_html(std::ostream* p_html, const std::string& dir, const std::string& filename_prefix, unsigned int exp_index) const;
  //! Load from file
  static Result load(std::istream* in);
};

void print_basic_morphism(std::ostream* p_out, const EndomorphismSLP<int>& e) {
  e.for_each_non_trivial_image([&p_out] (const EndomorphismSLP<int>::symbol_image_pair_type& pair) {
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


template<typename T>
void Result<T>::save(std::ostream* p_out) const {
  std::ostream& out = *p_out;

  out << EXPERIMENT_DELIMITER << std::endl;

  write_comment(p_out, "rank");
  out << rank_ << std::endl;

  write_comment(p_out, "original morphism value");
  out << values_.morphism_value << std::endl;

  write_comment(p_out, "conjugation value");
  out << values_.conjugation_value << std::endl;

  write_comment(p_out, "minimized morphism value");
  out << values_.minimized_morphism_value << std::endl;

  write_comment(p_out, "minimized conjugation value");
  out << values_.minimized_conjugation_value << std::endl;

  write_comment(p_out, "num of composed morphisms");
  out << composition_parts_.size() << std::endl;

  write_comment(p_out, "composed morphisms");
  for (const auto& morphism: composition_parts_) {
    morphism.save_to(p_out);
    write_comment(p_out, "");
  }

  write_comment(p_out, "num of conjugation morphisms");
  out << conjugation_parts_.size() << std::endl;

  write_comment(p_out, "conjugation morphisms");
  for (const auto& morphism: conjugation_parts_) {
    morphism.save_to(p_out);
    write_comment(p_out, "");
  }

  write_comment(p_out, "original morphism");
  morphism_.save_to(p_out);

  write_comment(p_out, "conjugation");
  conjugation_.save_to(p_out);

  write_comment(p_out, "minimized morphism");
  values_.minimized_morphism.save_to(p_out);

  write_comment(p_out, "minimized conjugation");
  values_.minimized_conjugation.save_to(p_out);

  write_comment(p_out, "num of morphism minimizing conjugators");
  out << values_.morphism_minimizing_sequence.size() << std::endl;

  write_comment(p_out, "sequence minimizing morphism");
  for (const auto& morphism: values_.morphism_minimizing_sequence) {
    morphism.save_to(p_out);
    write_comment(p_out, "");
  }

  write_comment(p_out, "num of conjugation minimizing conjugators");
  out << values_.conjugation_minimizing_sequence.size() << std::endl;

  write_comment(p_out, "sequence minimizing conjugaion");
  for (const auto& morphism: values_.conjugation_minimizing_sequence) {
    morphism.save_to(p_out);
    write_comment(p_out, "");
  }
}

template<typename T>
void Result<T>::print(std::ostream* p_out) const {
  std::ostream& out = *p_out;

  out << EXPERIMENT_DELIMITER << std::endl;

  write_comment(p_out, "rank");
  out << rank_ << std::endl;

  write_comment(p_out, "original morphism value");
  out << values_.morphism_value << std::endl;

  write_comment(p_out, "conjugation value");
  out << values_.conjugation_value << std::endl;

  write_comment(p_out, "minimized morphism value");
  out << values_.minimized_morphism_value << std::endl;

  write_comment(p_out, "minimized conjugation value");
  out << values_.minimized_conjugation_value << std::endl;

  write_comment(p_out, "num of composed morphisms");
  out << composition_parts_.size() << std::endl;

  write_comment(p_out, "composed morphisms");
  for (const auto& morphism: composition_parts_) {
    print_basic_morphism(p_out, morphism);
  }

  write_comment(p_out, "num of conjugation morphisms");
  out << conjugation_parts_.size() << std::endl;

  write_comment(p_out, "conjugation morphisms");
  for (const auto& morphism: conjugation_parts_) {
    print_basic_morphism(p_out, morphism);
  }

  write_comment(p_out, "original morphism");
  morphism_.save_graphviz(p_out, "Morphism");

  write_comment(p_out, "conjugation");
  conjugation_.save_graphviz(p_out, "Conjugation");

  write_comment(p_out, "minimized morphism");
  values_.minimized_morphism.save_graphviz(p_out, "MinMorphism");

  write_comment(p_out, "minimized conjugation");
  values_.minimized_conjugation.save_graphviz(p_out, "MinConjugation");

  write_comment(p_out, "num of morphism minimizing conjugators");
  out << values_.morphism_minimizing_sequence.size() << std::endl;

  write_comment(p_out, "sequence minimizing morphism");
  for (const auto& morphism: values_.morphism_minimizing_sequence) {
    print_basic_morphism(p_out, morphism);
  }

  write_comment(p_out, "num of conjugation minimizing conjugators");
  out << values_.conjugation_minimizing_sequence.size() << std::endl;

  write_comment(p_out, "sequence minimizing conjugaion");
  for (const auto& morphism: values_.conjugation_minimizing_sequence) {
    print_basic_morphism(p_out, morphism);
  }
}

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


template<typename T>
void Result<T>::print_html(std::ostream* p_html, const std::string& dir, const std::string& aux_filename_prefix, unsigned int exp_index) const {
  std::ostream& html = *p_html;

  html << indent(1) << "<h3>" << "Experiment with rank = " << rank_ << "</h3>" << std::endl;

  html << indent(1) << "<p>" << "Values:" << "</p>" << std::endl;
  html << indent(1) << "<ul>" << std::endl;
  html << indent(2) << "<li>Original morphism: " << values_.morphism_value << "</li>" << std::endl;
  html << indent(2) << "<li>Conjugation: " << values_.conjugation_value << "</li>" << std::endl;
  html << indent(2) << "<li>Minimized morphism: " << values_.minimized_morphism_value << "</li>" << std::endl;
  html << indent(2) << "<li>Minimized conjugation: " << values_.minimized_conjugation_value << "</li>" << std::endl;
  html << indent(1) << "</ul>" << std::endl;

  html << indent(1) << "<p>" << "num of composed morphisms: " << composition_parts_.size() << "</p>" << std::endl;
  html << indent(1) << "<p>" << "num of conjugation morphisms: " << conjugation_parts_.size() << "</p>" << std::endl;
  html << indent(1) << "<p>" << "num of conjugators minimizing the original morphism: " << values_.morphism_minimizing_sequence.size() << "</p>" << std::endl;
  html << indent(1) << "<p>" << "num of conjugators minimizaing conjugation: " << values_.conjugation_minimizing_sequence.size() << "</p>" << std::endl;



  html << indent(1) << "<table>" << std::endl;

  html << indent(2) << "<thead>" << std::endl;

  html << indent(3) << "<tr>" << std::endl;
  html << indent(4) << "<th class=\"left\" colspan=\"2\">Sample genrating data</th>" << std::endl;
  html << indent(4) << "<th colspan=\"2\">Minimizing conjugators</th>" << std::endl;
  html << indent(3) << "</tr>" << std::endl;

  html << indent(3) << "<tr>" << std::endl;
  html << indent(4) << "<th class=\"left\">Morphism</th>" << std::endl;
  html << indent(4) << "<th>Conjugators</th>" << std::endl;
  html << indent(4) << "<th>original</th>" << std::endl;
  html << indent(4) << "<th>conjgation</th>" << std::endl;
  html << indent(3) << "</tr>" << std::endl;

  html << indent(2) << "</thead>" << std::endl;

  html << indent(2) << "<tbody>" << std::endl;
  auto comp_iterator = composition_parts_.begin();
  auto conj_iterator = conjugation_parts_.begin();
  auto min_iterator = values_.morphism_minimizing_sequence.rbegin();
  auto min_conj_iterator = values_.conjugation_minimizing_sequence.rbegin();

  while (true) {
    bool all_done = true;
    html << indent(3) << "<tr>" << std::endl;
    if(comp_iterator != composition_parts_.end()) {
      all_done = false;
      html << indent(4) << "<td>";
      print_basic_morphism(p_html, *comp_iterator);
      html << "</td>" << std::endl;
      ++comp_iterator;
    } else {
      html << "<td></td>" << std::endl;
    }

    if(conj_iterator != conjugation_parts_.end()) {
      all_done = false;
      html << indent(4) << "<td>";
      print_basic_morphism(p_html, *conj_iterator);
      html << "</td>" << std::endl;
      ++conj_iterator;
    } else {
      html << "<td></td>" << std::endl;
    }

    if(min_iterator != values_.morphism_minimizing_sequence.rend()) {
      all_done = false;
      html << indent(4) << "<td>";
      print_basic_morphism(p_html, *min_iterator);
      html << "</td>" << std::endl;
      ++min_iterator;
    } else {
      html << "<td></td>" << std::endl;
    }

    if(min_conj_iterator != values_.conjugation_minimizing_sequence.rend()) {
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

  auto generate_morphism_description = [&] (const EndomorphismSLP<int>& e, const std::string& name, const std::string& description) {
    std::stringstream s;
    s << aux_filename_prefix << "_" << name << "_" << exp_index;
    std::string filename = s.str();
    std::string description_filename(filename + ".gv");
    std::string image_filename(filename + ".gif");

    html << indent(1) << "<p>" << description << "</p>";
    html << indent(1) << "<img src=\"" << image_filename << "\"/>" << std::endl;

    std::ofstream description_file(dir + description_filename);
    e.save_graphviz(&description_file, name);
    s.str("");
    s << "dot -Tgif " << dir << description_filename << " -o " << dir << image_filename;
    std::system(s.str().c_str());
  };

  generate_morphism_description(morphism_, "morphism", "Original morphism.");
  generate_morphism_description(conjugation_, "conjugation", "Conjugation.");
  generate_morphism_description(values_.minimized_morphism, "min_morphism", "Minimized morphism");
  generate_morphism_description(values_.minimized_conjugation, "min_conjugation", "Minimized conjugation.");
}

std::istream& operator>>(std::istream& in, LongInteger& n) {
  std::string s;
  std::getline(in, s);
  n = LongInteger(s);
}

template<typename T>
Result<T> Result<T>::load(std::istream* p_in) {
  std::istream& in = *p_in;

  Result result;
  skip_comments(p_in);
  in >> result.rank_;

  skip_comments(p_in);

//  std::string s;
//  std::getline(in, s);
//  std::cout << s;

  in >> result.values_.morphism_value;

  skip_comments(p_in);
  in >> result.values_.conjugation_value;

  skip_comments(p_in);
  in >> result.values_.minimized_morphism_value;

  skip_comments(p_in);
  in >> result.values_.minimized_conjugation_value;

  skip_comments(p_in);
  unsigned int composition_num;
  in >> composition_num;
  result.composition_parts_.reserve(composition_num);
  skip_comments(p_in);
  for (unsigned int i = 0; i < composition_num; ++i) {
    result.composition_parts_.push_back(EndomorphismSLP<int>::load_from(p_in));
    skip_comments(p_in);
  }

  unsigned int conjugation_num;
  in >> conjugation_num;
  result.conjugation_parts_.reserve(conjugation_num);
  for (unsigned int i = 0; i < conjugation_num; ++i) {
    skip_comments(p_in);
    result.conjugation_parts_.push_back(EndomorphismSLP<int>::load_from(p_in));
  }

  skip_comments(p_in);
  result.morphism_ = EndomorphismSLP<int>::load_from(p_in);

  skip_comments(p_in);
  result.conjugation_ = EndomorphismSLP<int>::load_from(p_in);

  skip_comments(p_in);
  result.values_.minimized_morphism = EndomorphismSLP<int>::load_from(p_in);

  skip_comments(p_in);
  result.values_.minimized_conjugation = EndomorphismSLP<int>::load_from(p_in);

  skip_comments(p_in);
  unsigned int min_morphism_conjugators_num;
  in >> min_morphism_conjugators_num;
  result.values_.morphism_minimizing_sequence.reserve(min_morphism_conjugators_num);
  skip_comments(p_in);
  for (unsigned int i = 0; i < min_morphism_conjugators_num; ++i) {
    result.values_.morphism_minimizing_sequence.push_back(EndomorphismSLP<int>::load_from(p_in));
    skip_comments(p_in);
  }

  skip_comments(p_in);
  unsigned int min_conjugation_conjugators_num;
  in >> min_conjugation_conjugators_num;
  result.values_.conjugation_minimizing_sequence.reserve(min_conjugation_conjugators_num);
  skip_comments(p_in);
  for (unsigned int i = 0; i < min_conjugation_conjugators_num; ++i) {
    result.values_.conjugation_minimizing_sequence.push_back(EndomorphismSLP<int>::load_from(p_in));
    skip_comments(p_in);
  }

  return result;
}


template <typename T>
std::ostream& operator<<(std::ostream& out, const TargetValues<T>& values) {
  return out << values.morphism_value << ";" << values.conjugation_value << ";"
                << values.minimized_morphism_value << ";" << values.minimized_conjugation_value << ";";
}



void conjugation_length_based_attack_statistics() {
  std::ofstream out("lba_based_on_total_length_result.txt");
  write_comment(&out, "Legnth-base attack to Conjugation Search Problem for Automorphisms of Free Group");
  write_comment(&out, "minimization value name");
  out << "total image legnth" << std::endl;

  typedef unsigned int uint;
  for (auto rank : {3}) {
    std::cout << "rank=" << rank << std::endl;
    UniformAutomorphismSLPGenerator<int> rnd(rank);
    for (auto size : {5}) {
      std::cout << "num of composed automorphisms=" << size << std::endl;
      for (auto conj_num: {size}) {
        std::cout << "num of conjugators=" << conj_num << std::endl;
        const uint iterations_num = 1;

        auto start_time = our_clock::now();
        for (int i = 0; i < iterations_num; ++i) {
          std::cout << "Iteration " << i << std::endl;

          auto iteration_start_time = our_clock::now();

          Result<LongInteger> result(rank, size, conj_num, rnd);
          result.morphism_ = minimize_morphism(rank, result.morphism_, &slp_vertices_num<int>).first;

          result.values_ = conjugation_function_minimization_attack(rank, result.morphism_, result.conjugation_,
            &slp_vertices_num<int>,//&total_images_length,//set the function which value we minimize when use conjugation attack
            true);
          result.save(&out);
          auto iteration_time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(our_clock::now() - iteration_start_time);
          write_comment(&out, "time");
          out << COMMENT_LINE_START << iteration_time_in_ms.count() <<  "ms" << std::endl;
        }
        auto time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(our_clock::now() - start_time);
        std::cout << "(iterations num=" << iterations_num << ",size=" << size << ",conjug num=" << conj_num << ",rank=" << rank << "): "
                  << time_in_ms.count() << "ms, "
                  << time_in_ms.count() / iterations_num << "ms per iteration, "
                  << std::endl;
      }
    }
  }
}


void lba_success_precentage() {
  std::ofstream out("lba_success_precentage_2.txt");
  write_comment(&out, "Legnth-base attack to Conjugation Search Problem for Automorphisms of Free Group");
  write_comment(&out, "minimization value name");
  out << "total image legnth" << std::endl;


  auto target_function = total_images_length;

  typedef unsigned int uint;
  for (auto rank : {3}) {
    std::cout << "rank=" << rank << std::endl;
    UniformAutomorphismSLPGenerator<int> rnd(rank);
    for (uint size: {2}) {
      std::cout << "num of composed automorphisms=" << size << std::endl;
      for (auto conj_num: {size}) {
        std::cout << "num of conjugators=" << conj_num << std::endl;
        const uint iterations_num = 3;

        auto start_time = our_clock::now();
        unsigned int success_num = 0;
        for (unsigned int i = 0; i < iterations_num; ++i) {
//          if (i % 10 == 0)
//            std::cout << "Iteration " << i << std::endl;

          auto iteration_start_time = our_clock::now();

          Result<LongInteger> result(rank, size, conj_num, rnd);
          result.morphism_ = minimize_morphism(rank, result.morphism_, target_function).first;

          result.values_ = conjugation_function_minimization_attack(rank, result.morphism_, result.conjugation_,
            target_function,//set the function which value we minimize when use conjugation attack
            false);
          if (result.values_.minimized_morphism_value == result.values_.minimized_conjugation_value &&
              result.values_.minimized_morphism == result.values_.minimized_conjugation)
            ++success_num;
          result.save(&out);
          auto iteration_time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(our_clock::now() - iteration_start_time);
          write_comment(&out, "time");
          out << COMMENT_LINE_START << iteration_time_in_ms.count() <<  "ms" << std::endl;
        }
        auto time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(our_clock::now() - start_time);
        std::cout << "(iterations num=" << iterations_num << ",size=" << size << ",conjug num=" << conj_num << ",rank=" << rank << "): "
                  << time_in_ms.count() << "ms, "
                  << time_in_ms.count() / iterations_num << "ms per iteration, "
                  << "success num=" << success_num
                  << std::endl;
      }
    }
  }
}

void skip_line(std::istream* in) {
  in->ignore(std::numeric_limits<std::streamsize>::max(),'\n');
}

template<typename T>
struct Results {
    std::string value_name;
    std::vector<Result<T> > exp_results;
};

template<typename T>
Results<T> read_results(const std::string& filename) {
  std::ifstream in(filename);
  in.exceptions (std::ios::eofbit | std::ios::failbit |
                 std::ios::badbit);
  skip_comments(&in);

  Results<T> results;
  std::getline(in, results.value_name);


  try {
  while (in && !in.eof())
      results.exp_results.push_back(Result<LongInteger>::load(&in));
  } catch(std::ifstream::failure e) {
    //finished reading
  }

  return results;
}

template<typename T>
void print_all_html(const std::string& dir, const std::string& filenames_prefix, const Results<T>& results) {
  print_html(dir, filenames_prefix, results, [] (const Result<T>& r) {return true;});
}

template<typename T, typename Filter>
void print_html(const std::string& dir, const std::string& filenames_prefix, const Results<T>& results, Filter filter = [] (const Result<T>& r) {return true;}) {
  std::string stylesheet_filename = filenames_prefix + "_stylesheet.css";

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


  std::ofstream html(dir + filenames_prefix + "_index.html");
  html << "<!DOCTYPE html>" << std::endl;
  html << "<html>" << std::endl;
  html << "<head>" << std::endl;
  html << "<link type=\"text/css\" rel=\"stylesheet\" href=\"" << stylesheet_filename << "\" />" << std::endl;
  html << "  <title>" << "Results for minimization value " << results.value_name << "</title>" << std::endl;
  html << "</head>" << std::endl;

  html << "<body>" << std::endl;
  int n = 0;
  for (const Result<T>& result: results.exp_results) {
    if (filter(result)) {
      result.print_html(&html, dir, filenames_prefix, n++);
    }
  }
  html << "</body>" << std::endl;

  html << "</html>" << std::endl;

}


int main() {
  //  composition_statistics();
//  conjugation_length_based_attack_statistics();
  lba_success_precentage();
  std::string filename("lba_success_precentage_2.txt");
  auto result = read_results<LongInteger>(filename);
  print_all_html("result/", "result", result);
}
