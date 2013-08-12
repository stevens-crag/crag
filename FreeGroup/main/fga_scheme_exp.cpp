/*
* fga_scheme_exp.cpp
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
std::default_random_engine rnd;

struct Stats {
    unsigned int height;
    unsigned int vertices_num;
};

bool logging = true;

namespace crag {

namespace fga_crypto {

class FGAExperiment {
  public:
    FGAExperiment() : FGAExperiment(&std::cout) {
    }

    FGAExperiment(std::ostream* p_out) : out_(*p_out) {
      out_ << "|u|,|v|,|c|,time,height,vertices_num" << std::endl;
    }

    void evaluate_scheme_time(const SchemeParameters& params, unsigned int samples_num) {
      std::cout << "Starting fga experiment with params = ("
                << params.A_SIZE << ", "
                << params.B_SIZE << ", "
                << params.U_LENGTH << ", "
                << params.V_LENGTH << ", "
                << params.C_SIZE << ")" << std::endl;
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
        print(s.height);
        print(s.vertices_num, false);

        std::cout << time_in_ms.count() << "ms" << std::endl;
      }
      auto s = ms / 1000;
      std::cout << "total time " << s << "s for " << samples_num << " samples (average time="
                << (s / samples_num) << "s)" <<  std::endl;
    }

    Stats evaluate_sample(const SchemeParameters& params) {
      KeysGenerator alice(params, rnd);
      KeysGenerator bob(params, rnd);

      alice.is_logging = true;

      auto a_pk = alice.public_keys();
      auto b_pk = bob.public_keys();

      auto a_processed = bob.process_incoming_public_keys(a_pk);

      auto a_shared_key = alice.make_shared_key(a_processed, true);
      auto key = a_shared_key();

      Stats s;
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

}// namespace fga_crypto
}// namespace crag


int main(int argc, char* argv[]) {
  if (argc != 2) {
    std::cout << "Wrong input arguments" << std::endl;
  } else {
    std::string output_filename(argv[1]);
    std::cout << "Writing to " << output_filename << std::endl;
    std::ofstream out(output_filename);
    crag::fga_crypto::FGAExperiment e(&out);
    int iterations = 5;
    for (int size = 5; size < 20; ++size) {
      e.evaluate_scheme_time(crag::fga_crypto::SchemeParameters(3, size, size, size), iterations);
    }
  }
}
