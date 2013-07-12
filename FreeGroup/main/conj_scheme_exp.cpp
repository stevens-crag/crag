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
#include "AAGCrypto.h"

typedef std::chrono::high_resolution_clock our_clock;
using namespace crag;
using namespace crag::fga_crypto;

//random engine for this instance
std::default_random_engine rnd;

struct Stats {
    LongInteger total_length;
    unsigned int height;
    unsigned int vertices_num;
};

namespace crag {

namespace fga_crypto {

class FGAExperiment {
  public:
    FGAExperiment() : FGAExperiment(&std::cout) {
    }

    FGAExperiment(std::ostream* p_out) : out_(*p_out) {
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

}// namespace fga_crypto

namespace aag_crypto {

class AAGExperiment {

  public:
    AAGExperiment() : AAGExperiment(&std::cout) {
    }

    AAGExperiment(std::ostream* p_out) : out_(*p_out) {
      out_ << "key_length,time,total_length,height,vertices_num," << std::endl;
    }

    void evaluate_scheme_time(const SchemeParameters& params, unsigned int samples_num) {
      unsigned long ms = 0;
      for (int i = 0; i < samples_num; ++i) {
        auto start_time = our_clock::now();
        Stats s = evaluate_sample(params);
        auto time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(our_clock::now() - start_time);
        ms += time_in_ms.count();

        print(params.KEY_LENGTH);
        print(time_in_ms.count());
        print(s.total_length);
        print(s.height);
        print(s.vertices_num);

        std::cout << time_in_ms.count() << "ms" << std::endl;
      }
      std::cout << "total time " << ms << "ms for " << samples_num << " samples" << std::endl;
    }

    Stats evaluate_sample(const SchemeParameters& params) {
      auto k_gen = make_keys_generator(params, &rnd);

      auto a_pub = k_gen.generate_public_key(Alice);
      auto b_pub = k_gen.generate_public_key(Bob);

      auto a_priv = k_gen.generate_private_key(a_pub);
      auto b_priv = k_gen.generate_private_key(b_pub);

      TransmittedInfo b_ti(a_pub, b_priv);

      auto key = k_gen.make_shared_key(a_priv, b_ti, Alice);

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

}// namespace aag_crypto

}// namespace crag


int main() {
  std::ofstream out("conj_exp.csv");
  crag::aag_crypto::AAGExperiment e(&out);
  int iterations = 10;
  for (int i = 0; i < 100; ++i) {
    e.evaluate_scheme_time(crag::aag_crypto::SchemeParameters(3, 20, 20, 5, 8, i), iterations);
  }
}
