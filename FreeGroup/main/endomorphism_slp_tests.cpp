
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



Statistic<LongInteger> get_endomorphism_images_lengths_stat(const EndomorphismSLP<int>& e) {
  Statistic<LongInteger> length_stat;
  auto length_calculator = [&] (const EndomorphismSLP<int>::symbol_image_pair_type& pair) {
    slp::Vertex v = pair.second;
    length_stat(v.length());
  };
  e.for_each_non_trivial_image(length_calculator);
  return length_stat;
}

void print_stats(std::ostream* out, const EndomorphismSLP<int>& e) {
  *out << "height=" << height(e) << ", vertices num=" << slp_vertices_num(e)
       << ", image lengths=(" << get_endomorphism_images_lengths_stat(e) << ")";
}

LongInteger total_images_length(const EndomorphismSLP<int>& em) {
  return sum(get_endomorphism_images_lengths_stat(em));
}

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

template<typename Func>
EndomorphismSLP<int> minimize_morphism(unsigned int rank, const EndomorphismSLP<int>& e, Func f, bool logging = true) {
  auto start_morphism = e.free_reduction();
  auto start_value = f(start_morphism);
  if (logging) {
    std::cout << "start value = " << start_value << ", ";
    print_stats(&std::cout, start_morphism);
    std::cout << std::endl;
  }
  while(true) {
    auto min_value = start_value;
    EndomorphismSLP<int> min_trial = EndomorphismSLP<int>::identity();

    EndomorphismSLP<int>::for_each_basic_morphism(rank,
                                       [&] (const EndomorphismSLP<int>& e) {
      auto trial = e * start_morphism * e.inverse();
      auto reduced_trial = trial.free_reduction();
      auto value = f(reduced_trial);
      if (value < min_value) {
        min_value = value;
        min_trial = reduced_trial;
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
    start_morphism = min_trial;
    start_value = min_value;
  }
  return start_morphism;
}

template<typename T>
struct TargetValues {
  T morphism_value;
  T conjugation_value;
  T minimized_morphism_value;
  T minimized_conjugation_value;

  TargetValues(T morphism_val, T conjug_val, T min_morphism_val, T min_conjug_val)
    : morphism_value(morphism_val), conjugation_value(conjug_val),
      minimized_morphism_value(min_morphism_val), minimized_conjugation_value(min_conjug_val) {}

  static std::string legend() {
    return "morphism value;conjugation value;minimized morphism value;minimized conjugator value;";
  }
};

template <typename T>
std::ostream& operator<<(std::ostream& out, const TargetValues<T>& values) {
  return out << values.morphism_value << ";" << values.conjugation_value << ";"
                << values.minimized_morphism_value << ";" << values.minimized_conjugation_value << ";";
}

template<typename TargetFunc, typename TargetValueType = LongInteger>
TargetValues<TargetValueType> conjugation_length_based_attack(unsigned int rank, const EndomorphismSLP<int>& e,
                                     const EndomorphismSLP<int>& e_conjugated,
                                     TargetFunc target_func,
                                     bool detailed_logging = true) {
  auto reduced_e = e.free_reduction();
  auto reduced_conj = e_conjugated.free_reduction();
  if (detailed_logging)
    std::cout << "Minimizing conjugation..." << std::endl;

  auto min = minimize_morphism(rank, reduced_conj, target_func, detailed_logging);

  if (detailed_logging)
    std::cout << "Minimizing e itself..." << std::endl;
  auto min_e = minimize_morphism(rank, reduced_e, target_func, detailed_logging);

  if (detailed_logging)
    std::cout << "Comparing the conjugation, e, and minimized versions." << std::endl;

  auto reduced_min = min.free_reduction();
  auto reduced_min_e = min_e.free_reduction();

  auto e_val = target_func(reduced_e);
  auto conj_val = target_func(reduced_conj);
  auto min_e_val = target_func(min_e);
  auto min_val = target_func(min);

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
  print_stats(&std::cout, reduced_min);
  std::cout << ", val=" << min_val << std::endl;



//  slp::MatchingTable mt;
  if (detailed_logging) {
    std::cout << "Images of reduced morphisms" << std::endl;
    for (int i = 1; i <= rank; ++i) {
      std::cout << "image of " << i << std::endl;
      auto img_e = reduced_min_e.image_word(i);
      auto img_min = reduced_min.image_word(i);
      std::cout << "e min :" << img_e << std::endl;
      std::cout << "min   :" << img_min << std::endl;
    }
  }

  return TargetValues<TargetValueType>(e_val, conj_val, min_e_val, min_val);
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

void conjugation_length_based_attack_statistics() {
  std::ofstream out("lba_based_on_total_length_result.csv");
  //printing head
  out << "rank;automorphism size;conjugators num;"
            << TargetValues<LongInteger>::legend()
            << std::endl;
  typedef unsigned int uint;
  for (auto rank : {3}) {
    std::cout << "rank=" << rank << std::endl;
    UniformAutomorphismSLPGenerator<int> rnd(rank);
    for (auto size : {10}) {
      std::cout << "num of composed automorphisms=" << size << std::endl;
      for (auto conj_num: {size}) {
        std::cout << "num of conjugators=" << conj_num << std::endl;
        const uint iterations_num = 1;

        auto start_time = our_clock::now();
        for (int i = 0; i < iterations_num; ++i) {
          std::cout << "Iteration " << i << std::endl;
          EndomorphismSLP<int> e = EndomorphismSLP<int>::composition(size, rnd);

          auto e_conjugation = e.conjugate_with(conj_num, rnd);

          auto result = conjugation_length_based_attack(rank, e, e_conjugation,
            &total_images_length,
            false);
          out << rank << ";" << size << ";" << conj_num << ";"
                    << result << std::endl;
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



int main() {
  //  composition_statistics();
  conjugation_length_based_attack_statistics();
}
