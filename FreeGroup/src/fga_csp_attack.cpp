/*
 * fga_csp_attack.cpp
 *
 *  Created on: July 19, 2013
 *      Author: pmorar
 */

#include "fga_csp_attack.h"


namespace crag {
namespace fga_csp_attack {

Statistic<LongInteger> get_endomorphism_images_lengths_stat(const Aut& e) {
  Statistic<LongInteger> length_stat;
  std::map<Aut::TerminalSymbol, LongInteger> lengths = images_length(e);
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
  std::map<Aut::TerminalSymbol, LongInteger> lengths = images_length(aut);
  for (const auto& pair: lengths) {
    sum += pair.second * pair.second;
  }
  return sqrt(sum);
}

LongInteger sqrt_of_sq_dif_of_lengths(const Aut& aut, const Aut& aut1) {
  LongInteger sum(0);
  std::map<Aut::TerminalSymbol, LongInteger> lengths = images_length(aut);
  std::map<Aut::TerminalSymbol, LongInteger> lengths1 = images_length(aut1);
  std::set<Aut::TerminalSymbol> parsed_keys;
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

const char COMMENT_LINE_START = '#';
const char* EXPERIMENT_DELIMITER = "#-------------------------";
const char* MORPHISMS_DELIMITER = "***";

void skip_comments(std::istream* in) {
  while (in->peek() == COMMENT_LINE_START ||
         in->peek() == '\n') {
    in->ignore(std::numeric_limits<std::streamsize>::max(),'\n');
  }
}

void skip_line(std::istream* in) {
  in->ignore(std::numeric_limits<std::streamsize>::max(),'\n');
}

//! Prints automorphism parameters.
void print_stats(std::ostream* out, const Aut& e) {
  *out << "height=" << height(e) << ", vertices num=" << slp_vertices_num(e)
       << ", image lengths=(" << get_endomorphism_images_lengths_stat(e) << ")";
}

//! Print stats for the interval [begin, end).
template<typename In>
void print_stats(std::ostream* out, In begin, In end) {
  *out << "{" << std::endl;
  std::for_each(begin, end, [&] (const Aut& aut) {
   print_stats(out, aut);
   *out << std::endl;
  });
  *out << "}" << std::endl;
}


void print_basic_morphism(std::ostream* p_out, const Aut& e) {
  e.for_each_non_trivial_image([&p_out] (const Aut::symbol_image_pair_type& pair) {
    *p_out << pair.first << " -> ";
    slp::Vertex v = pair.second;
    if (v.height() <= 1)
      *p_out << slp::TerminalVertex(v).terminal_symbol();
    else {
      auto l = v.left_child();
      auto r = v.right_child();
      *p_out << slp::TerminalVertex(l).terminal_symbol() << " "
                << slp::TerminalVertex(r).terminal_symbol();
    }
    *p_out << std::endl;
  });
}

void print_explicit_images(std::ostream* p_out, const Aut& aut) {
  std::ostream& out = *p_out;
  aut.for_each_non_trivial_image([&out] (const Aut::symbol_image_pair_type& pair) {
    out << "<div>" << std::endl;
    out << pair.first << " -> ";
    slp::VertexWord image(pair.second);
    std::for_each(image.begin(), image.end(), [&out] (int symbol) {out << symbol;});
    out << std::endl;
    out << "</div>" << std::endl;
  });
}


MinimizationResult::MinimizationResult(std::istream* p_in) {
  std::istream& in = *p_in;
  skip_comments(p_in);
  in >> min_value_;
  skip_comments(p_in);
  int num;
  in >> num;
  min_morphisms_.reserve(num);
  for (unsigned int i = 0; i < num; ++i) {
    skip_comments(p_in);
    min_morphisms_.push_back(Aut::load_from(p_in));
  }
  skip_comments(p_in);
  in >> num;

  minimizing_sequence_.reserve(num);
  for (unsigned int i = 0; i < num; ++i) {
    skip_comments(p_in);
    minimizing_sequence_.push_back(Aut::load_from(p_in));
  }
}

void MinimizationResult::save(std::ostream* p_out) const {
  std::ostream& out = *p_out;
  write_comment(p_out, "min value");
  out << min_value_ << std::endl;
  write_comment(p_out, "minimized elements");
  write_comment(p_out, "num");
  out << min_morphisms_.size() << std::endl;
  write_comment(p_out, "");
  for (const auto& morphism: min_morphisms_) {
    morphism.save_to(p_out);
    write_comment(p_out, "");
  }

  write_comment(p_out, "minimizing sequence");
  write_comment(p_out, "num");
  out << minimizing_sequence_.size() << std::endl;
  write_comment(p_out, "");
  for (const auto& morphism: minimizing_sequence_) {
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



void Result::print_html(std::ostream* p_html, const std::string& dir, const std::string& aux_filename_prefix, unsigned int exp_index, GenerateImages is_generating_images) const {
  std::ostream& html = *p_html;

  html << indent(1) << "<h3>" << "Experiment " << exp_index << "</h3>" << std::endl;

  html << indent(1) << "<p>" << "Parameters:" << "</p>" << std::endl;
  html << indent(1) << "<ul>" << std::endl;
  html << indent(2) << "<li>tuple size =" << aut_num_ << "</li>" << std::endl;
  html << indent(2) << "<li>rank = " << rank_ << "</li>" << std::endl;
  html << indent(2) << "<li>|base| = " << comp_num_ << "</li>" << std::endl;
  html << indent(2) << "<li>|conjugator| =" << conj_num_ << "</li>" << std::endl;
  html << indent(2) << "<li>|minimizing conjugator| = " << min_conjugations_.minimizing_sequence_.size() << "</li>" << std::endl;
  html << indent(1) << "</ul>" << std::endl;

  html << indent(1) << "<p>" << "Values:" << "</p>" << std::endl;
  html << indent(1) << "<ul>" << std::endl;
  html << indent(2) << "<li>morphism (minimized) = " << minimized_morphisms_.min_value_ << "</li>" << std::endl;
  html << indent(2) << "<li>conjugation = " << conjugation_value_ << "</li>" << std::endl;
  html << indent(2) << "<li>minimized conjugation = " << min_conjugations_.min_value_ << "</li>" << std::endl;
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
  auto min_conj_iterator = min_conjugations_.minimizing_sequence_.cbegin();

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

    if(min_conj_iterator != min_conjugations_.minimizing_sequence_.cend()) {
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

  Aut minimizator = Aut::composition(min_conjugations_.minimizing_sequence_.cbegin(), min_conjugations_.minimizing_sequence_.cend());
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
    morphisms_.push_back(AutDecomposition(p_in));
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


void print_values_to_csv_file(const string& filename, const Results& results) {
  std::ofstream csv_out(filename);
  csv_out << "val;conj_val;min_val" << std::endl;

  std::for_each(results.exp_results.cbegin(), results.exp_results.cend(), [&csv_out] (const Result& result) {
    csv_out << result.minimized_morphisms_.min_value_ << ";" <<
      result.conjugation_value_ << ";" <<
      result.min_conjugations_.min_value_ << std::endl;
  });
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



void KickResult::print_html(std::ostream* p_html) const {
  std::ostream& html = *p_html;

  html << indent(1) << "<h3>" << "Kick attempt</h3>" << std::endl;

  html << indent(1) << "<p>" << "Parameters:" << "</p>" << std::endl;
  html << indent(1) << "<ul>" << std::endl;
  html << indent(2) << "<li>tuple size =" << aut_num_ << "</li>" << std::endl;
  html << indent(2) << "<li>|conjugator| =" << conj_num_ << "</li>" << std::endl;
  html << indent(2) << "<li>|minimizing conjugator| = " << min_conjugations_.minimizing_sequence_.size() << "</li>" << std::endl;
  html << indent(1) << "</ul>" << std::endl;

  html << indent(1) << "<p>minimized value = " << min_conjugations_.min_value_ << "</p>" << std::endl;

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
  auto min_conj_iterator = min_conjugations_.minimizing_sequence_.cbegin();

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

    if(min_conj_iterator != min_conjugations_.minimizing_sequence_.cend()) {
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

void print_all_html(const std::string& dir, const std::string& filenames_prefix, const Results& results, GenerateImages is_generating_images) {
  print_html(dir, filenames_prefix, results, is_generating_images, [] (const Result& r) {return true;});
}

void print_not_successful_html(const std::string& dir, const std::string& filenames_prefix, const Results& results, GenerateImages is_generating_images) {
  print_html(dir, filenames_prefix, results, is_generating_images, [] (const Result& r) {
    return r.minimized_morphisms_.min_value_ != r.min_conjugations_.min_value_;//TODO equalities
  });
}




} //namespace fga_csp_attack
} //namespace crag
