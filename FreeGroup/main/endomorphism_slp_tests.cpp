/*
 * endomorphism_slp_tests.cpp
 *
 *  Created on: Mar 19, 2013
 *      Author: pmorar
 */

#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <cstdlib>

#include <sys/stat.h>
#include <pwd.h>

#include "EndomorphismSLP.h"

using namespace crag;

typedef std::chrono::high_resolution_clock our_clock;
typedef EndomorphismSLP<int> Aut;

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

int main(int argc, char* argv[]) {
  auto myuid = getuid();
  auto mypasswd = getpwuid(myuid);
  std::string dir(mypasswd->pw_dir);
  dir += "/Documents/exp_results/";
  normal_form_statistics(dir + "normal_form_stat_large.csv");
}
