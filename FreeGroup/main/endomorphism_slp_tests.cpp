
/*
 * endomorphism_slp_tests.cpp
 *
 *  Created on: Mar 19, 2013
 *      Author: pmorar
 */
#include <iostream>
#include <chrono>
#include "EndomorphismSLP.h"

typedef std::chrono::high_resolution_clock our_clock;
using namespace crag;

template<typename T>
class SimpleStat {
  public:
    static_assert(std::is_integral<T>::value, "SimpleStat parameter must be of integral type");

    SimpleStat()
      : sum_(0),
        min_(std::numeric_limits<T>::max()),
        max_(std::numeric_limits<T>::min()),
        count_(0) {}

    void reset() {
      sum_ = 0;
      min_ = std::numeric_limits<T>::max();
      max = std::numeric_limits<T>::min();
      count_ = 0;
    }

    SimpleStat& operator+=(const SimpleStat& s) {
      sum_ += s.sum_;
      if (min_ > s.min_)
        min_ = s.min_;
      if (max_ < s.max_)
        max_ = s.max_;
      count_ += s.count_;
    }

    void add_value(T val) {
      sum_ += val;
      if (val < min_)
        min_ = val;
      if (val > max_)
        max_ = val;
      ++count_;
    }

    long double average() const {
      return static_cast<long double>(sum_) / count_;
    }

    T max() const {
      return max_;
    }

    T min() const {
      return min_;
    }

  private:
    long long int sum_;
    T min_;
    T max_;
    unsigned long long int count_;
};

template<>
class SimpleStat<LongInteger> {
  public:
    typedef LongInteger Int;

    SimpleStat()
      : sum_(0),
        min_(0),
        max_(0),
        count_(0) {}

    void reset() {
      sum_ = 0;
      min_ = 0;
      max_ = 0;
      count_ = 0;
    }

    SimpleStat& operator+=(const SimpleStat& s) {
      sum_ += s.sum_;
      if (count_ == zero()) {
        min_ = s.min_;
        max_ = s.max_;
      } else {
        if (min_ > s.min_)
          min_ = s.min_;
        if (max_ < s.max_)
          max_ = s.max_;
      }
      count_ += s.count_;
    }

    void add_value(const Int& val) {
      sum_ += val;
      if (count_ > zero()) {
        if (val < min_)
          min_ = val;
        if (val > max_)
          max_ = val;
      } else {
        min_ = max_ = val;
      }
      ++count_;
    }

    Int average() const {
      return sum_ / count_;
    }

    Int max() const {
      return max_;
    }

    Int min() const {
      return min_;
    }

  private:
    static Int zero() {
      static Int zero_ = Int();
      return zero_;
    }

    Int sum_;
    Int min_;
    Int max_;
    Int count_;
};

template <typename T>
std::ostream& operator<<(std::ostream& out, const SimpleStat<T>& stat) {
  return out << "(avg = " << stat.average() << ", "
                 << "min = " << stat.min() << ", "
                 << "max = " << stat.max() << ")";
}

void composition_statistics() {
  typedef unsigned int uint;
  std::cout << "Legend: (num iterations, num of composed elements, rank)" << std::endl;
  for (auto rank : {3, 5, 10, 20}) {
    UniformAutomorphismSLPGenerator<int> rnd(rank);
    for (auto size : {100, 1000, 2000}) {
      const uint iterations_num = 100;
      SimpleStat<unsigned int> heightStat;
      SimpleStat<unsigned int> verticesNumStat;

      our_clock::duration time;
      for (int i = 0; i < iterations_num; ++i) {
        auto start_time = our_clock::now();
        auto e = EndomorphismSLP<int>::composition(size, rnd);
        time += our_clock::now() - start_time;
        auto height = crag::height(e);
        auto vertices_num = slp_vertices_num(e);

        heightStat.add_value(height);
        verticesNumStat.add_value(vertices_num);
      }
      auto time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(time);
      std::cout << "(iterations=" << iterations_num << ", size=" << size << ", rank=" << rank << "): "
                << time_in_ms.count() << "ms, "
                << time_in_ms.count() / iterations_num << "ms per iteration, "
                << "height " << heightStat << ", "
                << "vertices num " << verticesNumStat
                << std::endl;
    }
  }
}

SimpleStat<LongInteger> get_endomorphism_images_lengths_stat(const EndomorphismSLP<int>& e) {
  SimpleStat<LongInteger> length_stat;
  auto length_calculator = [&] (const EndomorphismSLP<int>::symbol_image_pair_type& pair) {
    slp::Vertex v = pair.second;
    length_stat.add_value(v.length());
  };
  e.for_each_non_trivial_image(length_calculator);
  return length_stat;
}

void print_stats(std::ostream* out, const EndomorphismSLP<int>& e) {
  *out << "height=" << height(e) << ", vertices num=" << slp_vertices_num(e)
       << ", image lengths=(" << get_endomorphism_images_lengths_stat(e) << ")";
}

template<typename Func>
void enumerate_elementary_automorphisms(int rank, Func f) {
  for (int i = 1; i <= rank; ++i)
    f(EndomorphismSLP<int>::inverter(i));
  for (int i = 1; i <= rank; ++i)
    for (int j = -rank; j <= rank; ++j)
      if (j != i && j != -i && j != 0) {
        f(EndomorphismSLP<int>::left_multiplier(j, i));
        f(EndomorphismSLP<int>::right_multiplier(i, j));
      }
}

SimpleStat<LongInteger> image_length_distance(unsigned int rank, const EndomorphismSLP<int>& e1,
                                              const EndomorphismSLP<int>& e2) {
  SimpleStat<LongInteger> length_stat;
  for (int i = 1; i <= rank; ++i) {
    auto img1 = e1.image(i);
    auto img2 = e2.image(i);
    LongInteger d = img1.length() - img2.length();
    length_stat.add_value(d > 0 ? d : -d);
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

    enumerate_elementary_automorphisms(rank,
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

void conjugation_length_based_attack(unsigned int rank, const EndomorphismSLP<int>& e,
                                     const EndomorphismSLP<int>& e_conjugated, bool detailed_logging = true) {
  auto reduced_conj = e_conjugated.free_reduction();
  if (detailed_logging)
    std::cout << "Minimizing conjugation..." << std::endl;

  auto target_func = [&] (const EndomorphismSLP<int>& em) -> LongInteger {
    return get_endomorphism_images_lengths_stat(em).max();
//    return image_length_distance(rank, em, e).average();
  };

  auto min = minimize_morphism(rank, reduced_conj, target_func, detailed_logging);

  if (detailed_logging)
    std::cout << "Minimizing e itself..." << std::endl;
  auto min_e = minimize_morphism(rank, e, target_func, detailed_logging);

  if (detailed_logging)
    std::cout << "Comparing the conjugation, e, and minimized versions." << std::endl;

  auto reduced_min = min.free_reduction();
  auto reduced_min_e = min_e.free_reduction();
  std::cout << "e :        ";
  print_stats(&std::cout, e);
  std::cout << std::endl;

  std::cout << "conj :     ";
  print_stats(&std::cout, reduced_conj);
  std::cout << std::endl;

  std::cout << "min e:     ";
  print_stats(&std::cout, reduced_min_e);
  std::cout << std::endl;

  std::cout << "min conj : ";
  print_stats(&std::cout, reduced_min);
  std::cout << std::endl;



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
}


void conjugation_length_based_attack_statistics() {
  typedef unsigned int uint;
  for (auto rank : {3, 5}) {
    std::cout << "rank=" << rank << std::endl;
    UniformAutomorphismSLPGenerator<int> rnd(rank);
    for (auto size : {10, 15}) {
      std::cout << "num of composed automorphisms=" << size << std::endl;
      for (auto conj_num: {10, 15}) {
        std::cout << "num of conjugators=" << conj_num << std::endl;
        const uint iterations_num = 5;

        auto start_time = our_clock::now();
        for (int i = 0; i < iterations_num; ++i) {
          std::cout << "Iteration " << i << std::endl;
          EndomorphismSLP<int> e = EndomorphismSLP<int>::composition(size, rnd);

          auto e_conjugation = e.conjugate_with(conj_num, rnd);

          conjugation_length_based_attack(rank, e, e_conjugation, false);

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
