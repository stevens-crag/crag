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
typedef EndomorphismSLP Aut;

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
    UniformAutomorphismSLPGenerator<> rnd(rank);
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

template <class TPermutation>
class AlternativeBasePermutations {
  public:
    static const std::vector<TPermutation>& permutations() {
      static std::vector<TPermutation> permutations = {
        TPermutation(), //for null terminal
        TPermutation({1, 0, 3, 11, 12, 6, 7, 13, 9, 8, 15, 2, 4, 5, 10, 14}), //permutation of the maximal order in S16
        TPermutation({8, 3, 9, 15, 7, 0, 1, 12, 5, 13, 10, 14, 2, 6, 4, 11}), //combined with the previous, it can give the whole group
      };
      return permutations;
    }
};

typedef crag::slp::TVertexHashAlgorithms<
  crag::slp::hashers::ImageLengthHash,
  crag::slp::hashers::SinglePowerHash,
  crag::slp::hashers::PermutationHashBase<crag::Permutation16, AlternativeBasePermutations<crag::Permutation16>>
> AlternativeVertexHashAlgorithms;


//! Difference between the vertices number and the number of hashes of vertices (hash given by the template parameter).
template<typename VertexHashAlgorithms = endomorphism_default_parameters::WeakVertexHashAlgorithms>
std::size_t v_num_h_num_gap(const EndomorphismSLP& e) {
  typename VertexHashAlgorithms::Cache cache;
  std::unordered_set<slp::Vertex> visited_vertices;
  std::unordered_set<typename VertexHashAlgorithms::VertexHash> visited_hashes;

  auto acceptor = [&visited_vertices] (const slp::inspector::InspectorTask& task) {
    return visited_vertices.count(task.vertex) == 0
        && visited_vertices.count(task.vertex.negate()) == 0;
  };

  auto inspect_root =[&acceptor,&visited_vertices,&visited_hashes,&cache] (const typename EndomorphismSLP::symbol_image_pair_type& v) {
    slp::Inspector<slp::inspector::Postorder, decltype(acceptor)> inspector(v.second, acceptor);
    while (!inspector.stopped()) {
      auto v = inspector.vertex();
      visited_vertices.insert(v);
      visited_hashes.insert(VertexHashAlgorithms::get_vertex_hash(v, &cache));
      inspector.next();
    }
  };

  e.for_each_non_trivial_image(inspect_root);
  return visited_vertices.size() - visited_hashes.size();
}


void hash_collisions_statistics(std::ostream* p_out, const uint rank, const uint size, const uint iterations) {
  std::ostream& out = *p_out;
  UniformAutomorphismSLPGenerator<> rnd(rank);
  uint total_collisions_num = 0;
  uint collisions_num = 0;
  uint gap_num = 0;
  uint total_gap = 0;
  auto time = our_clock::now();
  for (uint i = 0; i < iterations; ++i) {
    out << i << ":";
    out.flush();
    const auto e = Aut::composition(size, rnd);
    out << "gen,";
    out.flush();
    const auto regular = v_num_h_num_gap(e);

    std::size_t alt = 0;
    try {
      out << "h1,";
      out.flush();
      if (regular > 0) {
        out << "pos gap, ";
        out.flush();
        ++gap_num;
        total_gap += regular;
      }
      alt = v_num_h_num_gap<AlternativeVertexHashAlgorithms>(e);
      out << "h2,";
      out.flush();

    } catch(const std::exception& exception) {
      std::cerr << "Exception!"  << std::endl;
      throw exception;
    } catch(...) {
      std::cerr << "Unidentified error!" << std::endl;
      throw std::exception();
    }

    auto n = regular - alt;
    if (n != 0) {
      n = n > 0 ? n : -n;
      total_collisions_num += n;
      ++collisions_num;
      e.save_to(p_out);
      out << " Collision! ";
      out.flush();
    }
    out << "done. ";
    if (i % 4 == 3)
      out << std::endl;
    out.flush();
  }
  auto time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(our_clock::now() - time);
  out << std::endl << "Summary: rank=" << rank
      << ",size=" << size
      << ",iterations=" << iterations
      << ",collisions=" << collisions_num
      << ",total collisions=" << total_collisions_num
      << ",gap num=" << gap_num
      << ",total gap=" << total_gap
      << ",time=" << time_in_ms.count() << std::endl;
}

const std::string TAB("    ");
const std::string MODE("--mode=");
const std::string NF("nf");
const std::string COLLISIONS("collisions");

const std::string RANK("--rank=");
const std::string SIZE("--size=");
const std::string ITERATIONS("--iter=");

void print_usage() {
  std::cout << "Options:" << std::endl;
  std::cout << MODE << " : mode of calculations." << std::endl;
  std::cout << TAB << "Possible values:" << std::endl;
  std::cout << TAB << NF << " : calculating normal form statistics." << std::endl;
  std::cout << TAB << COLLISIONS << " : calculating hash collisions statistics." << std::endl;
  std::cout << RANK << " : automorphisms rank." << std::endl;
  std::cout << SIZE << " : automorphisms composition size." << std::endl;
  std::cout << ITERATIONS << " : iterations num." << std::endl;
}




int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cout << "Not enough parameters." << std::endl;
    print_usage();
  }
  std::string mode;
  uint rank = 0;
  uint size = 0;
  uint iterations = 0;
  for (int i = 1; i < argc; ++i) {
    std::string option(argv[i]);
    if (option.find(MODE) == 0) {
      mode = option.substr(MODE.size());
    } else if (option.find(RANK) == 0) {
      rank = std::stoi(option.substr(RANK.size()));
    } else if (option.find(SIZE) == 0) {
      size = std::stoi(option.substr(SIZE.size()));
    } else if (option.find(ITERATIONS) == 0) {
      iterations = std::stoi(option.substr(ITERATIONS.size()));
    }
  }

  if (mode == NF) {
    auto myuid = getuid();
    auto mypasswd = getpwuid(myuid);
    std::string dir(mypasswd->pw_dir);
    dir += "/Documents/exp_results/";
    //TODO refactor normal_form_statistics to accept
    normal_form_statistics(dir + "normal_form_stat_large.csv");
  } else if (mode == COLLISIONS) {
    if(rank == 0 || size == 0 || iterations == 0) {
      print_usage();
      return 0;
    }
    hash_collisions_statistics(&std::cout, rank, size, iterations);
  } else {
    print_usage();
  }

}
