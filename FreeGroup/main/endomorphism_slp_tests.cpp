
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

    void addValue(T val) {
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

        heightStat.addValue(height);
        verticesNumStat.addValue(vertices_num);
      }
      auto time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(time);
      std::cout << "(num=" << iterations_num << ", size=" << size << ", rank=" << rank << "): "
                << time_in_ms.count() << "ms, "
                << time_in_ms.count() / iterations_num << "ms per iteration, "
                << "height " << heightStat << ", "
                << "vertices num " << verticesNumStat
                << std::endl;
    }
  }
}


int main() {
  composition_statistics();
}
