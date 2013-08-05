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

bool logging = true;

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
      out_ << "key_length,time,height,vertices_num," << std::endl;
    }

    void evaluate_scheme_time(const SchemeParameters& params, CalculationType calc_type, unsigned int samples_num) {
      if (logging) {
        std::cout << "Starting with params = (" << params.ALICE_TUPLE_SIZE << ", "
                  << params.BOB_TUPLE_SIZE << ", "
                  << params.LOWER_BOUND_PUB_KEY_LENGTH << ", "
                  << params.UPPER_BOUND_PUB_KEY_LENGTH << ", "
                  << params.KEY_LENGTH << ")" << std::endl;
      }
      unsigned long ms = 0;
      for (int i = 0; i < samples_num; ++i) {
        auto start_time = our_clock::now();
        Stats s = evaluate_sample(params, calc_type);
        auto time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(our_clock::now() - start_time);
        ms += time_in_ms.count();

        print(params.KEY_LENGTH);
        print(time_in_ms.count());
        print(s.height);
        print(s.vertices_num, false);

        std::cout << time_in_ms.count() << "ms" << std::endl;
      }
      std::cout << "total time " << ms << "ms for " << samples_num << " samples" << std::endl;
    }

    Stats evaluate_sample(const SchemeParameters& params, CalculationType calc_type) {
      auto k_gen = make_keys_generator(params, &rnd);

      auto a_pub = k_gen.generate_public_key(Alice);
      auto b_pub = k_gen.generate_public_key(Bob);

      if (logging) {
        std::cout << "pub key sizes" << std::endl << "A = (";
        for (int i = 0; i < params.ALICE_TUPLE_SIZE; ++i) {
          std::cout << slp_vertices_num(a_pub[i]()) << ", ";
        }
        std::cout << ")" << std::endl;

        std::cout << "B = (";
        for (int i = 0; i < params.BOB_TUPLE_SIZE; ++i) {
          std::cout << slp_vertices_num(b_pub[i]()) << ", ";
        }
        std::cout << ")" << std::endl;
      }

      auto a_priv = k_gen.generate_private_key(a_pub);
      auto b_priv = k_gen.generate_private_key(b_pub);

      if (logging) {
        std::cout << "private key sizes A = " << slp_vertices_num(a_priv()()) <<
                     ", B = " << slp_vertices_num(b_priv()()) << std::endl;
      }

      TransmittedInfo b_ti(a_pub, b_priv);

      if (logging) {
        std::cout << "transmitted info sizes" << std::endl;
        std::cout << "B = (";
        for (int i = 0; i < params.BOB_TUPLE_SIZE; ++i) {
          std::cout << slp_vertices_num(b_ti[i]()) << ", ";
        }
        std::cout << ")" << std::endl;
      }
      Aut key = k_gen.make_shared_key(a_priv, b_ti, Alice, calc_type);

      if (logging) {
        std::cout << "shared key size " << slp_vertices_num(key) << std::endl;
      }

      Stats s;
      s.total_length = 0;
//      auto lengths = images_length(key);
//      for (const auto& pair: lengths) {
//        s.total_length += pair.second;
//      }
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
        out_ << "," << std::endl;
      }
    }
};

}// namespace aag_crypto

}// namespace crag


int main(int argc, char* argv[]) {
  if (argc != 3) {
    std::cout << "Wrong input arguments" << std::endl;
  } else {
    std::string calc_type(argv[1]);
    crag::aag_crypto::CalculationType type;
    if (calc_type == "--block") {
      type = crag::aag_crypto::BlockFrNf;
      std::cout << "Block mode" << std::endl;
    } else if (calc_type == "--iterative") {
      type = crag::aag_crypto::IterativeFrNf;  
      std::cout << "Iterative mode" << std::endl;
    } else {
      type = crag::aag_crypto::SinlgeFrNf;  
      std::cout << "Single mode" << std::endl;
    }

    std::string output_filename(argv[2]);
    std::cout << "Writing to " << output_filename << std::endl;
    std::ofstream out(output_filename);
    crag::aag_crypto::AAGExperiment e(&out);
    int iterations = 5;
    for (int i = 10; i < 100; i+=10) {
      e.evaluate_scheme_time(crag::aag_crypto::SchemeParameters(3, 20, 20, 4, 5, i), type, iterations);
    }
  }
}