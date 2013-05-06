/*
 * FGACrypto.h
 *
 *  Created on: Apr 27, 2013
 *      Author: pmorar
 */

#ifndef CRAG_FREE_GROUPS_FGACRYPTO_H
#define CRAG_FREE_GROUPS_FGACRYPTO_H

#include "EndomorphismSLP.h"

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

      SchemeParameters(int n, int a, int b, int u, int v, int c)
        : N(n), A_SIZE(a), B_SIZE(b), U_LENGTH(u), V_LENGTH(v), C_SIZE(c) {}

      const unsigned int N;
      const unsigned int A_SIZE;
      const unsigned int B_SIZE;
      const unsigned int U_LENGTH;
      const unsigned int V_LENGTH;
      const unsigned int C_SIZE;
      //TODO make constructor
  };

  typedef AutomorphismDescription<EndomorphismSLP<int> > AutomorphismDescription;

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

      CommutatorSet conjugate_with(const AutomorphismDescription& conjugator) const{
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


      AutomorphismDescription get_s(int index) const {
        if (index > 0) {
          return s[index - 1];
        } else {
          return s[- index - 1].inverse_description();
        }
      }

      AutomorphismDescription get_r(int first_index, int second_index) const {
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

//        for (int i = 0; i < 4; ++i) {
//          std::cout << "alpha " << i + 1 << std::endl;
//          alphas_[i]().print(&std::cout);
//        }

//        for (int i = 0; i < 4; ++i) {
//          std::cout << "beta " << i + 1 << std::endl;
//          betas_[i]().print(&std::cout);
//        }

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

//        std::cout << "u ";
//        for (int i: u_)
//          std::cout << i << " ";
//        std::cout << std::endl;

//        std::cout << "v ";
//        for (int i: v_)
//          std::cout << i << " ";
//        std::cout << std::endl;

        //generating public and private keys bases
        s_.reserve(4);
        for (int i = 0; i < 4; ++i) {
          s_.push_back(alphas_[i] * betas_[i]);
        }

//        for (int i = 0; i < 4; ++i) {
//          std::cout << "s " << i + 1 << std::endl;
//          s_[i]().print(&std::cout);
//        }

        r_.reserve(4);
        for (int i: {0, 1})
          for (int j: {2, 3}) {
//            std::cout << "r_" << i + 1 << "_" << j + 1 << std::endl;
            r_.push_back(CommutatorSet(betas_[i], betas_[j]));
//            r_.back().get(false, false)().print(&std::cout);
          }

        priv_key_base_ = make_commutator(get_betas_composition(u_),
            get_betas_composition(v_));

//        std::cout << "priv key base" << std::endl;
//        priv_key_base_().print(&std::cout);

        //generating keys
        priv_key_ = priv_key_base_.conjugate_with(c_);

//        std::cout << "priv key" << std::endl;
//        priv_key_().print(&std::cout);

        pub_keys_.s = conjugate_all(s_, c_);
        pub_keys_.r = conjugate_all(r_, c_);
      }

      const PublicKeys& public_keys() const {
        return pub_keys_;
      }

      const AutomorphismDescription& private_key() const {
        return priv_key_;
      }

      //! Processes public keys provided by other party to send them back
      PublicKeys process_incoming_public_keys(const PublicKeys& incoming_public_keys) {
        PublicKeys result;
        result.s = conjugate_all(incoming_public_keys.s, priv_key_);
        result.r = conjugate_all(incoming_public_keys.r, priv_key_);
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
        for (int row_index: v_) {
          key *= calculate_private_key_line(row_index, processed_public_keys).conjugate_with(conjugator);
          conjugator *= processed_public_keys.get_s(row_index);//TODO optimize
        }

        if (order) {
          key = priv_key_ * key.inverse_description();//Alice: a * (bab^-1)^-1
        } else {
          key *= priv_key_.inverse_description();//Bob: aba^-1 *= b^-1
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
        //TODO cache
        AutomorphismDescription conjugator;
        AutomorphismDescription value;
        for (int col_index: u_) {
          auto beta_conj = public_key.get_r(col_index, row_index);
          value = beta_conj.conjugate_with(conjugator) * value;
          conjugator *= public_key.get_s(col_index);
        }
        return value;
      }

   };
} // namespace fga_crypto
} // namespace crag

#endif // CRAG_FREE_GROUPS_FGACRYPTO_H
