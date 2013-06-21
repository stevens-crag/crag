/*
 * FGACrypto.h
 *
 *  Created on: Apr 27, 2013
 *      Author: pmorar
 */

#ifndef CRAG_FREE_GROUPS_FGACRYPTO_H
#define CRAG_FREE_GROUPS_FGACRYPTO_H

#include "EndomorphismSLP.h"

#include <fstream>

/**
 * The file contains cryptoscheme programmed corresponding to the specifications by Sasha Ushakov.
 */

namespace crag {

namespace fga_crypto {

  struct SchemeParameters {
      static const SchemeParameters canonical() {
        static SchemeParameters params = SchemeParameters(3, 5, 4, 25, 25, 50);
//        N = 3;
//        A_SIZE = 5;
//        B_SIZE = 4;
//        U_LENGTH = 25;
//        V_LENGTH = 25;
//        C_SIZE = 50;
        return params;
      }


      SchemeParameters(int n, int u, int v, int c)
        : N(n), A_SIZE(5), B_SIZE(4), U_LENGTH(u), V_LENGTH(v), C_SIZE(c) {}

      SchemeParameters(int n, int a, int b, int u, int v, int c)
        : N(n), A_SIZE(a), B_SIZE(b), U_LENGTH(u), V_LENGTH(v), C_SIZE(c) {}


      const unsigned int N;
      const unsigned int A_SIZE;
      const unsigned int B_SIZE;
      const unsigned int U_LENGTH;
      const unsigned int V_LENGTH;
      const unsigned int C_SIZE;
  };

  typedef AutomorphismDescription<EndomorphismSLP<int> > AutomorphismDescription;



  void print_stats(const AutomorphismDescription& aut_d) {
    auto a = aut_d();
    std::cout << "vertices num (a=" << slp_vertices_num(a);
    auto fr = a.free_reduction();
    std::cout << ", fr=" << slp_vertices_num(fr);
    auto nf = fr.normal_form();
    std::cout << ", nf=" << slp_vertices_num(nf) << ")" << std::endl;
  }

  //! Set of commutators for a given pair of automorphisms. We use it only because we can not find inverses efficiently.
  class CommutatorSet {
    public:
      CommutatorSet(const AutomorphismDescription& first, const AutomorphismDescription& second) {
        auto inv_first = first.inverse_description();
        auto inv_second = second.inverse_description();
        comm_.reserve(4);
        comm_.push_back(make_commutator(first, second));
        comm_.push_back(make_commutator(inv_first, second));
        comm_.push_back(make_commutator(first, inv_second));
        comm_.push_back(make_commutator(inv_first, inv_second));
      }

      const AutomorphismDescription& get(bool first_inversed, bool second_inversed) const {
        return comm_[(first_inversed ? 1 : 0) + (second_inversed ? 2 : 0)];
      }

      CommutatorSet conjugate_with(const AutomorphismDescription& conjugator) const {
        CommutatorSet result;
        result.comm_ = conjugate_all(comm_, conjugator);
        return result;
      }

    private:
      std::vector<AutomorphismDescription> comm_;

      CommutatorSet() {}
  };


  //! Public keys
  struct PublicKeys {
      std::vector<AutomorphismDescription> s;
      std::vector<CommutatorSet> r;

      PublicKeys() {}

      PublicKeys(std::vector<AutomorphismDescription>&& v_s, std::vector<CommutatorSet>&& v_r)
        : s(v_s), r(v_r) {}

      AutomorphismDescription get_s(int index) const {
        if (index > 0) {
          return s[index - 1];
        } else {
          return s[- index - 1].inverse_description();
        }
      }

      const AutomorphismDescription& get_r(int first_index, int second_index) const {
        bool first_inversed = first_index < 0;
        bool second_inversed = second_index < 0;
        first_index = first_inversed ? - first_index : first_index;
        first_index -= 1;
        second_index = second_inversed ? - second_index : second_index;
        second_index -= 3;
        return r[2 * first_index + second_index].get(first_inversed, second_inversed);
      }
  };

  //! Generates keys for the scheme.
  class KeysGenerator {
      static const bool debug = false;
    public:

      KeysGenerator() = delete;

      //! Create random keys according to the specified parameters and random generator.
      template<typename RandomEngine>
      KeysGenerator(const SchemeParameters& params, RandomEngine& rand)
        : params_(params),
          alphas_(),
          betas_(),
          u_(),
          v_(),
          s_(),
          r_(),
          priv_key_base_(),
          c_(),
          priv_key_(),
          pub_keys_(),
          shared_key_() {

        //generating alphas, betas
        UniformAutomorphismSLPGenerator<int, RandomEngine> random_for_alphas(params.N, &rand);
        UniformAutomorphismSLPGenerator<int, RandomEngine> random_for_betas(params.N + 1, 2 * params.N, &rand);
        UniformAutomorphismSLPGenerator<int, RandomEngine> random_for_c(2 * params.N, &rand);

        alphas_.reserve(4);
        betas_.reserve(4);
        for (int i = 0; i < 4; ++i) {
          alphas_.push_back(AutomorphismDescription(params.A_SIZE, random_for_alphas));
          betas_.push_back(AutomorphismDescription(params.B_SIZE, random_for_betas));
        }

        if (debug) {
          for (int i = 0; i < 4; ++i) {
            std::cout << "alpha " << i + 1 << " n=" << alphas_[i].composed_num() << std::endl;
            alphas_[i]().print(&std::cout);
          }

          for (int i = 0; i < 4; ++i) {
            std::cout << "beta " << i + 1 << " n=" << betas_[i].composed_num() <<  std::endl;
            betas_[i]().print(&std::cout);
          }
        }

        c_ = AutomorphismDescription(params.C_SIZE, random_for_c);

        //generating u, v
        std::uniform_int_distribution<int> binary_rand(0, 1);
        auto next_gen = [&binary_rand, &rand] (int shift) {
          int base = shift + binary_rand(rand);
          return binary_rand(rand) == 0 ? base : -base;
        };

        auto random_binary_vector = [&next_gen] (int size, int shift) {
          std::vector<int> v;
          v.reserve(size);
          int next_generator = next_gen(shift);
          v.push_back(next_generator);
          for (unsigned int i = 1; i < size; ++i) {
            int prev = v[i - 1];
            do {
              next_generator = next_gen(shift);
            } while (next_generator == -prev);
            v.push_back(next_generator);
          }
          return v;
        };

        u_ = random_binary_vector(params.U_LENGTH, 1);
        v_ = random_binary_vector(params.V_LENGTH, 3);

        if (debug) {
          std::cout << "u ";
          for (int i: u_)
            std::cout << i << " ";
          std::cout << std::endl;

          std::cout << "v ";
          for (int i: v_)
            std::cout << i << " ";
          std::cout << std::endl;
        }

        //generating public and private keys bases
        s_.reserve(4);
        for (int i = 0; i < 4; ++i) {
          s_.push_back(alphas_[i] * betas_[i]);
        }

        if (debug)
          for (int i = 0; i < 4; ++i) {
            std::cout << "s " << i + 1 << " n=" << s_[i].composed_num() <<  std::endl;
            s_[i]().print(&std::cout);
          }

        r_.reserve(4);
        for (int i: {0, 1})
          for (int j: {2, 3}) {
            r_.push_back(CommutatorSet(betas_[i], betas_[j]));
            if (debug) {
              std::cout << "r_" << i + 1 << "_" << j + 1 << " n=" << r_.back().get(false, false).composed_num() <<  std::endl;
              r_.back().get(false, false)().print(&std::cout);
            }
          }

        priv_key_base_ = make_commutator(get_betas_composition(u_),
            get_betas_composition(v_));

        if (debug) {
          std::cout << "priv key base" << " n=" << priv_key_base_.composed_num() <<   std::endl;
//          priv_key_base_().print(&std::cout);
          print_stats(priv_key_base_);
        }

        //generating keys
        priv_key_ = priv_key_base_.conjugate_with(c_);

        if (debug) {
          std::cout << "priv key" << " n=" << priv_key_.composed_num() <<  std::endl;
//          priv_key_().print(&std::cout);
          print_stats(priv_key_);
        }

        pub_keys_.s = conjugate_all(s_, c_);
        pub_keys_.r = conjugate_all(r_, c_);

//        std::cout << "pub key" << std::endl;
//        std::cout << "s" << std::endl;
//        for (auto& ad: pub_keys_.s) {
//          print_stats(ad);
//        }
//        std::cout << "r" << std::endl;
//        for (auto& cs: pub_keys_.r) {
//          print_stats(cs.get(false, false));
//        }
      }

      const PublicKeys& public_keys() const {
        return pub_keys_;
      }

      const AutomorphismDescription& private_key() const {
        return priv_key_;
      }

      //! Processes public keys provided by other party to send them back
      PublicKeys process_incoming_public_keys(const PublicKeys& incoming_public_keys) {
        PublicKeys result(conjugate_all(incoming_public_keys.s, priv_key_),
                            conjugate_all(incoming_public_keys.r, priv_key_));
//        std::cout << "processed pub key" << std::endl;
//        std::cout << "s" << std::endl;
//        for (auto& ad: result.s) {
//          print_stats(ad);
//        }
//        std::cout << "r" << std::endl;
//        for (auto& cs: result.r) {
//          print_stats(cs.get(false, false));
//        }
        return result;
      }

      //! Make shared keys with the given public keys of another party
      /**
       * @param processed_public_keys our public keys processed by another party
       * @param order if true Alice makes the key, otherwise Bob
       * @return
       */
      AutomorphismDescription make_shared_key(const PublicKeys& processed_public_keys, bool order = true) {
        AutomorphismDescription key;
        AutomorphismDescription conjugator;
        std::map<int, AutomorphismDescription> line_cache;
//        for (int row_index: v_) {
////          std::cout << "row " << row_index << std::endl;
//          auto cached_ad = line_cache.find(row_index);
//          AutomorphismDescription line;
//          if (cached_ad != line_cache.end()) {
//            line = cached_ad->second;
//          } else {
//            line = calculate_private_key_line(row_index, processed_public_keys);
//            line_cache.insert(std::make_pair(row_index, line));
//          }
//          key *= line.conjugate_with(conjugator);
//          conjugator *= processed_public_keys.get_s(row_index);
//        }


        for (int i = v_.size() - 1; i >= 0; --i) {
          const int row_index = v_[i];
          auto cached_ad = line_cache.find(row_index);
          AutomorphismDescription line;
          if (cached_ad != line_cache.end()) {
            line = cached_ad->second;
          } else {
            std::cout << "start line" << std::endl;
            line = calculate_private_key_line(row_index, processed_public_keys);
            std::cout << "line finished" << std::endl;
            line_cache.insert(std::make_pair(row_index, line));
          }
          key = line * key;
          if (i > 0) {
            const int conj_index = v_[i - 1];
            key = key.conjugate_with(processed_public_keys.get_s(conj_index));
          }
          std::cout << "key vert num (val=" << slp_vertices_num(key()) <<
                       ", inv=" << slp_vertices_num(key.inverse()) << ")" << std::endl;
          key = key.free_reduction().normal_form();
          std::cout << "nf key item vert num (val=" << slp_vertices_num(key()) <<
                       ", inv=" << slp_vertices_num(key.inverse()) << ")" << std::endl;
        }

        if (order) {
          key = priv_key_ * key.inverse_description();//Alice: a * (bab^-1)^-1
        } else {
          key *= priv_key_.inverse_description();//Bob: aba^-1 *= b^-1
        }
        std::cout << "almost final key vert num (val=" << slp_vertices_num(key()) <<
                     ", inv=" << slp_vertices_num(key.inverse()) << ")" << std::endl;
        key = key.free_reduction().normal_form();
        std::cout << "result key item vert num (val=" << slp_vertices_num(key()) <<
                     ", inv=" << slp_vertices_num(key.inverse()) << ")" << std::endl;

        if (debug) {
//          key().print(&std::cout);
        }
        return key;
      }

    private:
      const SchemeParameters params_;
      std::vector<AutomorphismDescription> alphas_;
      std::vector<AutomorphismDescription> betas_;
      std::vector<int> u_;
      std::vector<int> v_;
      std::vector<AutomorphismDescription> s_;
      std::vector<CommutatorSet> r_;

      AutomorphismDescription priv_key_base_;

      AutomorphismDescription c_;

      AutomorphismDescription priv_key_;

      PublicKeys pub_keys_;

      AutomorphismDescription shared_key_;

      AutomorphismDescription get_betas_composition(const std::vector<int>& pattern) {
        auto pick = [&] (int i) {
          if (i > 0) {
            return betas_[i - 1];
          } else {
            return betas_[- i - 1].inverse_description();
          }
        };
        AutomorphismDescription result;
        for (auto i: pattern) {
          result *= pick(i);
        }
        return result;
      }

      AutomorphismDescription calculate_private_key_line(int row_index, const PublicKeys& public_key) {
        AutomorphismDescription value;
//        for (int col_index: u_) {
////          std::cout << "column " << col_index << std::endl;
//          auto beta_conj = public_key.get_r(col_index, row_index);
//          value = beta_conj.conjugate_with(conjugator) * value;
//          conjugator *= public_key.get_s(col_index);
//        }

        for (int i = u_.size() - 1; i >= 0; --i) {
          const int col_index = u_[i];
          value *= public_key.get_r(col_index, row_index);
          const int conj_index = u_[i-1];
          if (i > 0) {
            value = value.conjugate_with(public_key.get_s(conj_index));
          }
          std::cout << "item before red vert num (val=" << slp_vertices_num(value()) <<
                      ", inv=" << slp_vertices_num(value.inverse()) << ")" << std::endl;
          value = value.free_reduction().normal_form();//reducing the size
          std::cout << "item vert num (val=" << slp_vertices_num(value()) <<
                       ", inv=" << slp_vertices_num(value.inverse()) << ")" << std::endl;
        }
        return value;
      }

   };
} // namespace fga_crypto
} // namespace crag

#endif // CRAG_FREE_GROUPS_FGACRYPTO_H
