/*
 * lba_attack.cpp
 *
 *  Created on: July 19, 2013
 *      Author: pmorar
 */


#include <cstdlib>

#include <sys/stat.h>
#include <pwd.h>


#include "fga_csp_attack.h"

using namespace crag;
using namespace crag::fga_csp_attack;

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


  MultiplicationsEnumerator enumerator(rank);
  UniformAutomorphismSLPGenerator<> rnd(rank);

  int eq_num = 0;
  auto start_time = our_clock::now();
  for (int i = 0; i < iterations_num; ++i) {
    std::cout << "iteration " << i << " ";

    auto iteration_start_time = our_clock::now();

    Result result(aut_num, rank, comp_num, conj_num, rnd, target_func, enumerator, &out);

    result.save(&out);

    auto iteration_time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(our_clock::now() - iteration_start_time);
    auto iteration_time_in_s = std::chrono::duration_cast<std::chrono::seconds>(our_clock::now() - iteration_start_time);
    write_comment(&out, "time");
    out << COMMENT_LINE_START << iteration_time_in_ms.count() <<  "ms" << std::endl;
    std::cout << iteration_time_in_s.count() << "s" << std::endl;
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



void process_length_base_attack() {
  //parameters
  typedef unsigned int uint;
  const uint aut_num = 3;
  const uint rank = 3;
  const uint comp_num = 10;
  const uint conj_num = 10;
  const uint iterations_num = 5;

  //filenames and dirs
//  auto myuid = getuid();
//  auto mypasswd = getpwuid(myuid);
//  std::string dir(mypasswd->pw_dir);
//  dir += "exp_res/";

  std::stringstream s;
  s << "lba_total_sum_a" << aut_num << "comp" << comp_num << "conj" << conj_num << "it" << iterations_num;
  std::string name(s.str());
  std::string filename = name + ".txt";
//  std::string csv_filename = dir + name + ".csv";
//  std::string html_dir = dir + name + "/";

  //work
  conjugation_length_based_attack_statistics(filename, aut_num, rank, comp_num, conj_num, iterations_num);
//  auto results = read_results(filename);

//  mode_t process_mask = umask(0);
//  int result_code = mkdir(html_dir.c_str(), S_IRWXU | S_IRWXG | S_IRWXO);
//  if (!result_code && errno == EEXIST) {
//    std::cerr << "can not create dir!" << std::endl;
//  }
//  umask(process_mask);

//  print_all_html(html_dir, "lba", results, not_generate_images);
//  print_values_to_csv_file(csv_filename, results);
}


void process_file_to_csv(const std::string& filename) {
  std::cout << "Reading results from " << filename << std::endl;
  auto results = read_results(filename);
  std::size_t separator_position = filename.find_last_of('.');
  auto csv_filename = filename.substr(0, separator_position) + ".csv";
  std::cout << "Writing results to " << csv_filename << "... ";
  print_values_to_csv_file(csv_filename, results);
  std::cout << "Done." << std::endl;
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
  auto min_conjs_ = result.min_conjugations_.min_morphisms_;

  AllNielsenGeneratorsEnumerator enumerator(rank);
  UniformAutomorphismSLPGenerator<> rnd(rank);

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


void kick_attempt() {
  const unsigned int conj_num = 10;
  const unsigned int iterations_num = 10;

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



void print_help() {
  std::cout << "Length based attack. Parameters:" << std::endl;
  std::cout << "--create_csv filename: create csv with summary from file." << std::endl;
}


int main(int argc, char* argv[]) {
  if (argc == 1) {
//    print_help();

    //  kick_attempt();
      process_length_base_attack();
    //  auto myuid = getuid();
    //  auto mypasswd = getpwuid(myuid);
    //  std::string dir(mypasswd->pw_dir);
    //  dir += "/Documents/exp_results/";

    //  normal_form_statistics(dir + "normal_form_stat_large.csv");
  } else {
    std::string command(argv[1]);
    if (command == "--create_csv") {
        std::string filename(argv[2]);
        process_file_to_csv(filename);
    } else {
        print_help();
    }

    /*
    std::string filename(argv[1]);
    process_file_to_csv(filename)*/;
  }
}
