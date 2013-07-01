/*
* conj_scheme_exp.cpp
*
* Created on: Jun 25, 2013
* Author: pmorar
*/
#include <iostream>
#include <fstream>
#include <sstream>
#include <chrono>

#include "FGACrypto.h"

typedef std::chrono::high_resolution_clock our_clock;
using namespace crag;
using namespace crag::fga_crypto;

//random engine for this instance
std::default_random_engine rnd(22345);


class Experiment {
    struct Stats {
        LongInteger total_length;
        unsigned int height;
        unsigned int vertices_num;
    };

  public:
    Experiment() : Experiment(&std::cout) {
    }

    Experiment(std::ostream* p_out) : out_(*p_out) {
      out_ << "|u|,|v|,|c|,time,total_length,height,vertices_num" << std::endl;
    }

    void evaluate_scheme_time(const SchemeParameters& params, unsigned int samples_num) {
      unsigned int n = 0;
      unsigned long ms = 0;
      for (int i = 0; i < samples_num; ++i) {
        auto start_time = our_clock::now();
        Stats s = evaluate_sample(params);
        auto time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(our_clock::now() - start_time);
        ms += time_in_ms.count();

        print(params.U_LENGTH);
        print(params.V_LENGTH);
        print(params.C_SIZE);
        print(time_in_ms.count());
        print(s.total_length);
        print(s.height);
        print(s.vertices_num, false);

        std::cout << time_in_ms.count() << "ms" << std::endl;
      }
      std::cout << "total time " << ms << "ms for " << samples_num << " samples" << std::endl;
    }

    Stats evaluate_sample(const SchemeParameters& params) {
      KeysGenerator alice(params, rnd);
      KeysGenerator bob(params, rnd);

      auto a_pk = alice.public_keys();
      auto b_pk = bob.public_keys();

      auto a_processed = bob.process_incoming_public_keys(a_pk);

      auto a_shared_key = alice.make_shared_key(a_processed, true);
      auto key = a_shared_key();

      Stats s;
      s.total_length = 0;
      auto lengths = images_length(key);
      for (const auto& pair: lengths) {
        s.total_length += pair.second;
      }
      s.height = height(key);
      s.vertices_num = slp_vertices_num(key);
      return s;
    }

  private:
    std::ostream& out_;

    template<typename T>
    void print(const T& val, bool separator = true) {
      out_ << val;
      if (separator) {
        out_ << ",";
      } else {
        out_ << std::endl;
      }
    }
};


int main() {
  int iterations = 10;
  std::ofstream out("conj_scheme3.csv");
  Experiment e(&out);
  { int base = 3;
    for (int i = 1; i <= 10; ++i) {
      e.evaluate_scheme_time(SchemeParameters(3, i, base, base), iterations);
      e.evaluate_scheme_time(SchemeParameters(3, base, i, base), iterations);
      e.evaluate_scheme_time(SchemeParameters(3, base, base, i), iterations);
    }
  }
//  crush_example();

}
