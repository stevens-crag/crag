
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
    LongInteger d = slp::reduce(img1).length() - slp::reduce(img2).length();
    length_stat.add_value(d > 0 ? d : -d);
  }
  return length_stat;
}

void conjugation_length_based_attack(unsigned int rank, const EndomorphismSLP<int>& e,
                                     const EndomorphismSLP<int>& e_conjugated) {
  SimpleStat<LongInteger> length_stat;
  auto length_calculator = [&] (const EndomorphismSLP<int>::symbol_image_pair_type& pair) {
    slp::Vertex v = pair.second;
    length_stat.add_value(slp::reduce(v).length());
  };

  length_stat.reset();
  e_conjugated.for_each_non_trivial_image(length_calculator);
  std::cout << "original conjugation: h=" << height(e_conjugated) << ", size=" << slp_vertices_num(e_conjugated)
            << ", l=(" << length_stat << ");" << std::endl;

  auto reduced = e_conjugated.free_reduction();
  length_stat.reset();
  reduced.for_each_non_trivial_image(length_calculator);
  std::cout << "reduced conjugation: h=" << height(reduced) << ", size=" << slp_vertices_num(reduced)
            << ", l=(" << length_stat << ");" << std::endl;

  std::cout << "Looking for the best deconjugators to minimize target function..." << std::endl;

  auto start_morphism = e_conjugated.free_reduction();
  auto start_value = image_length_distance(rank, e, start_morphism).max();
  std::cout << "initial value = " << start_value << std::endl;
  while(true) {
    auto min_value = start_value;
    EndomorphismSLP<int> min_trial = EndomorphismSLP<int>::identity();
    decltype(length_stat) min_stat;

    enumerate_elementary_automorphisms(rank,
                                       [&] (const EndomorphismSLP<int>& e) {
      auto trial = e * start_morphism * e.inverse();
      length_stat.reset();
      trial.for_each_non_trivial_image(length_calculator);
      auto value = image_length_distance(rank, e, trial).max();
      if (value < min_value) {
        min_value = value;
        min_trial = trial;
        min_stat = length_stat;
      }
    });

    if (start_value <= min_value) {
      std::cout << "Could not decrease target length." << std::endl;
      break;
    }
    std::cout << "Success: value=" << min_value << ", h=" << height(min_trial) << ", size=" << slp_vertices_num(min_trial)
              << ", l=(" << min_stat << ");" << std::endl;
    start_morphism = min_trial.free_reduction();
    start_value = min_value;
  }
  std::cout << "Comparing the minmized automorphism and e." << std::endl;

  length_stat.reset();
  e.for_each_non_trivial_image(length_calculator);
  std::cout << "e stats: h=" << height(e) << ", size=" << slp_vertices_num(e)
            << ", l=(" << length_stat << ");" << std::endl;

  length_stat.reset();
  start_morphism.for_each_non_trivial_image(length_calculator);
  std::cout << "result stats: h=" << height(start_morphism) << ", size=" << slp_vertices_num(start_morphism)
            << ", l=(" << length_stat << ");" << std::endl;

  auto free_reduced_result = start_morphism.free_reduction();
  length_stat.reset();
  free_reduced_result.for_each_non_trivial_image(length_calculator);
  std::cout << "free reduced result stats: h=" << height(free_reduced_result) << ", size=" << slp_vertices_num(free_reduced_result)
            << ", l=(" << length_stat << ");" << std::endl;

//  slp::MatchingTable mt;
  for (int i = 1; i <= rank; ++i) {
    std::cout << "image of " << i << std::endl;
    auto img_e = e.image_word(i);
    auto img_result = free_reduced_result.image_word(i);
    std::cout << "e   :" << img_e << std::endl;
    std::cout << "res :" << img_result << std::endl;

  }
}


void conjugation_length_based_attack_statistics() {
  typedef unsigned int uint;
  std::cout << "Legend: (num iterations, num of composed elements, rank, num_of_conjugators)" << std::endl;
//  std::vector<EndomorphismSLP<int>> inverses;
  for (auto rank : {3}) {
    UniformAutomorphismSLPGenerator<int> rnd(rank);
    for (auto size : {5}) {
//      inverses.reserve(size);
      for (auto conj_num: {1}) {
        const uint iterations_num = 1;

        auto start_time = our_clock::now();
        for (int i = 0; i < iterations_num; ++i) {
          std::cout << "Iteration " << i << std::endl;
          EndomorphismSLP<int> e = EndomorphismSLP<int>::composition(size, rnd);
//          inverses.clear();
//          for (int k = 0; k < size; ++k) {
//            auto next = rnd();
//            e *= next;
//            inverses.push_back(next.inverse());
//          }
//          EndomorphismSLP<int> e_inverse = EndomorphismSLP<int>::composition(inverses.rbegin(), inverses.rend());

          auto e_conjugation = e.conjugate_with(conj_num, rnd);

          conjugation_length_based_attack(rank, e, e_conjugation);

        }
        auto time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(our_clock::now() - start_time);
        std::cout << "(iterations=" << iterations_num << ",size=" << size << ",conjug=" << conj_num << ",rank=" << rank << "): "
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
