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
  public:
    Experiment() : Experiment(&std::cout) {
    }

    Experiment(std::ostream* p_out) : out_(*p_out) {
      out_ << "|u|,|v|,|c|,time,id_num" << std::endl;
    }

    void evaluate_scheme_time(const SchemeParameters& params, unsigned int samples_num) {
      unsigned int n = 0;
      unsigned long ms = 0;
      for (int i = 0; i < samples_num; ++i) {
        auto start_time = our_clock::now();
        n += evaluate_sample(params);
        auto time_in_ms = std::chrono::duration_cast<std::chrono::milliseconds>(our_clock::now() - start_time);
       ms += time_in_ms.count();
        std::cout << time_in_ms.count() << "ms" << std::endl;
      }
      std::cout << "total time " << ms << "ms for " << samples_num << " samples with " << n << " identities keys" << std::endl;
      std::vector<unsigned long int> a = {params.U_LENGTH, params.V_LENGTH, params.C_SIZE, ms, n};
      print_csv_line(a);
    }

    bool evaluate_sample(const SchemeParameters& params) {
      KeysGenerator alice(params, rnd);
      KeysGenerator bob(params, rnd);

      auto a_pk = alice.public_keys();
      auto b_pk = bob.public_keys();

      auto a_processed = bob.process_incoming_public_keys(a_pk);

      auto a_shared_key = alice.make_shared_key(a_processed, true);
      return a_shared_key().is_identity();
    }

  private:
    std::ostream& out_;

    template<typename T>
    void print_csv_line(const std::vector<T>& vals) {
      for (std::size_t i = 0; i < vals.size() - 1; ++i) {
        out_ << vals[i] << ",";
      }
      out_ << vals[vals.size() - 1] << std::endl;
    }
};

void crush_example() {
  std::ifstream in("crush_log");
  auto e = EndomorphismSLP<int>::load_from(&in);
  std::cout << "loaded" << std::endl;
  auto fr = e.free_reduction();
  std::cout << "fr num (val=" << slp_vertices_num(fr) << ")" << std::endl;
  std::cout << "free red" << std::endl;
  auto nf = fr.normal_form();
  std::cout << "normal form" << std::endl;
}


int main() {
  Experiment e;
  e.evaluate_scheme_time(SchemeParameters(3, 1, 1, 1), 5);
  e.evaluate_scheme_time(SchemeParameters(3, 2, 2, 2), 5);
//  crush_example();

}
