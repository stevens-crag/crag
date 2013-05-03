/*
 * FGACrypto.h
 *
 *  Created on: Jan 29, 2013
 *      Author: pmorar
 */

#ifndef CRAG_FREE_GROUPS_FGACRYPTO_H
#define CRAG_FREE_GROUPS_FGACRYPTO_H

#include "EndomorphismSLP.h"

/**
 * The file contains cryptoscheme programmed corresponding to the specifications by Sasha Ushakov.
 */

namespace crag {

  namespace FGACrypto {

    struct SchemeParameters {
        const unsigned int N = 3;
        const unsigned int A_SIZE = 5;
        const unsigned int B_SIZE = 4;
        const unsigned int U_LENGTH = 25;
        const unsigned int V_LENGTH = 25;
        const unsigned int C_SIZE = 50;
        //TODO make constructor
    };

    //! Set of commutators for a given pair of automorphisms. We use it only because we can not find inverses efficiently.
    class CommutatorSet {

        CommutatorSet(const AutomorphismDescription& first, const AutomorphismDescription& second) {
          comm_.reserve(4);
          comm_.push_back(AutomorphismDescription::make_commutator(first, second));
          comm_.push_back(AutomorphismDescription::make_commutator(first.inverse(), second));
          comm_.push_back(AutomorphismDescription::make_commutator(first, second.inverse()));
          comm_.push_back(AutomorphismDescription::make_commutator(first.inverse(), second.inverse()));
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
    };


    struct PublicKeys {
        std::vector<AutomorphismDescription> s;
        std::vector<CommutatorSet> r;


        AutomorphismDescription get_s(int index) {
          if (index > 0) {
            return s[index - 1];
          } else {
            return s[- index - 1].inverse_description();
          }
        }

        AutomorphismDescription get_r(int first_index, int second_index) {
          bool first_inversed = first_index < 0;
          bool second_inversed = second_index < 0;
          first_index = first_inversed ? - first_index : first_index;
          first_index -= 1;
          second_index = first_inversed ? - second_index : second_index;
          second_index -= 3;
          return r[2 * first_index + second_index].get(first_inversed, second_inversed);
        }
    };

    class KeysGenerator {
      public:

        KeysGenerator() = delete;

        //! Create random keys according to the specified parameters and random generator.
        template<typename RandomEngine>
        Keys(const SchemeParameters& params, RandomEngine& rand)
          : params_(params),
            alphas_(),
            betas_(),
            u_(),
            v_(),
            s_(),
            r_(),
            priv_key_base_(),
            c_(params.C_SIZE, UniformAutomorphismSLPGenerator<int, RandomEngine>(2 * params.N, p_rand)),
            priv_key_(),
            pub_keys_(),
            shared_key_() {

          //generating alphas, betas
          UniformAutomorphismSLPGenerator<int, RandomEngine> random_for_alphas(params.N, p_rand);
          UniformAutomorphismSLPGenerator<int, RandomEngine> random_for_betas(params.N + 1, 2 * params.N, p_rand);

          alphas_.reserve(4);
          betas_.reserve(4);
          for (int i = 0; i < 4; ++i) {
            alphas_.push_back(AutomorphismDescription(params.A_SIZE, random_for_alphas));
            betas_.push_back(AutomorphismDescription(params.B_SIZE, random_for_betas));
          }

          //generating u, v
          std::uniform_int_distribution<int> rand_generators(1, 2);
          u_.reserve(params.U_LENGTH);
          int next_generator = *p_rand() < 0.5 ? rand_generators(rand) : -rand_generators(rand);
          u_.push_back(next_generator);
          for (unsigned int i = 1; i < params.U_LENGTH; ++i) {
            int prev = u_[i - 1];
            do {
              next_generator = *p_rand() < 0.5 ? rand_generators(rand) : -rand_generators(rand);
            } while (next_generator != -prev);
            u_.push_back(next_generator);
          }

          v_.reserve(params.V_LENGTH);
          next_generator = *p_rand() < 0.5 ? params.N + rand_generators(rand) : -params.N - rand_generators(rand);
          v_.push_back(next_generator);
          for (unsigned int i = 1; i < params.V_LENGTH; ++i) {
            int prev = v_[i - 1];
            do {
              next_generator = *p_rand() < 0.5 ? params.N + rand_generators(rand) : -params.N - rand_generators(rand);
            } while (next_generator != -prev);
            v_.push_back(next_generator);
          }

          //generating public and private keys bases
          s_.reserve(4);
          for (int i = 0; i < 4; ++i) {
            s_.push_back(alphas_[i] * betas_[i]);
          }

          auto commutator = [] (const AutomorphismDescription& a1, const AutomorphismDescription& a2) {
            return (a1 * a2 * a1.inverse_description() * a2.inverse_description()).free_reduction();
          };

          r_.reserve(4);
          for (int i: {0, 1})
            for (int j: {3, 4})
              r_.push_back(CommutatorSet(betas_[i], betas_[j]));

          priv_key_base_ = commutator(get_betas_composition(u_),
              get_betas_composition(v_));

          //generating keys
          priv_key_ = a.conjugate_with(c);

          pub_keys_.s = conjugate_all(s_, c);
          pub_keys_.r = conjugate_all(r_, c);
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
        AutomorphismDescription make_shared_key(const PulbicKeys& processed_public_keys, bool order = true) {
          AutomorphismDescription key;
          AutomorphismDescription conjugator;
          for (int row_index: v_) {
            conjugator *= public_key.get_s(row_index);
            key *= calculate_private_key_line(row_index, processed_public_keys).conjugate_with(conjugator);
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

        static AutomorphismDescription get_betas_composition(const std::vector<int>& pattern) {
          auto pick = [&] (int i) {
            if (i > 0) {
              return betas_[i - 1]();
            } else {
              return betas_[- i - 1].inverse();
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
            conjugator *= public_key.get_s(col_index);
            auto beta_conj = public_key.get_r(col_index, row_index);
            value = beta_conj.conjugate_with(conjugator) * value;
          }
          return value;
        }

     };
  }
}

#endif // CRAG_FREE_GROUPS_FGACRYPTO_H
