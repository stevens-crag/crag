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

    //! A pair of automorphism and its inverse. We use it only because we can not find inverses effeciently.
    struct AutomorphismInvPair {
        EndomorphismSLP<int> a;
        EndomorphismSLP<int> a_inv;

        static AutomorphismInvPair make_conjugation(const AutomorphismDescription& morphism, const AutomorphismDescription& conjugator) {
          AutomorphismInvPair result;
          result.a = (conjugator.inverse() * morphism() * conjugator()).free_reduction();
          result.a_inv = (conjugator.inverse() * morphism.inverse() * conjugator()).free_reduction();
          return result;
        }

        static std::vector<AutomorphismInvPair> make_conjugations(const std::vector<AutomorphismDescription>& morphisms, const AutomorphismDescription& conjugator) {
          std::vector<AutomorphismInvPair> result;
          result.reserve(morphisms.size());
          for (int i = 0; i < morphisms.size(); ++i) {
            result.push_back(make_conjugation(morphisms[i], conjugator));
          }
          return result;
        }
    };

    struct PublicKeys {
        std::vector<AutomorphismInvPair> s;
        std::vector<AutomorphismInvPair> r;
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
          int next_generator = *p_rand() < 0.5 ? rand_generators(rand) : -rand_generators(rand);
          u_.push_back(next_generator);
          for (unsigned int i = 1; i < params.U_LENGTH; ++i) {
            int prev = u_[i - 1];
            do {
              next_generator = *p_rand() < 0.5 ? rand_generators(rand) : -rand_generators(rand);
            } while (next_generator != -prev);
            u_.push_back(next_generator);
          }

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
            return a1 * a2 * a1.invert() * a2.invert();
          };

          r_.reserve(4);
          for (unsigned int i = 0; i <= 1; ++i)
            for (unsigned int j = 3; j <= 4; ++j)
              r_.push_back(commutator(betas_[i], betas_[j]));

          priv_key_base_ = commutator(get_automorphism_composition(betas_[0], betas_[1], u_),
              get_automorphism_composition(betas_[2], betas_[3], v_));

          //generating keys
          priv_key_ = AutomorphismInvPair::make_conjugation(a, c);

          pub_keys_.s = AutomorphismInvPair::make_conjugations(s_, c);
          pub_keys_.r = AutomorphismInvPair::make_conjugations(r_, c);
        }

        const PublicKeys& public_keys() const {
          return pub_keys_;
        }

        const AutomorphismDescription& private_key() const {
          return priv_key_;
        }

        //! Processes public keys provided by other party to send them back
        PublicKeys process_incoming_public_keys(const PublicKeys& other_public_keys) {
          PublicKeys result;
          result.s = AutomorphismInvPair::make_conjugations(other_public_keys.s, c);
          result.r = AutomorphismInvPair::make_conjugations(other_public_keys.r, c);
          return result;
        }

        //! Make shared keys with the given public keys of another party
        /**
         * @param processed_public_keys our public keys processed by another party
         * @param order chooses the order in which we multiply our private key with the key made from processed public keys.
         * @return
         */
        AutomorphismInvPair make_shared_key(const PulbicKeys& processed_public_keys, bool order = true) {
          EndomorphismSLP<int> key;

          if (order) {
            key = priv_key_.a_inv * key;
          } else {
            key = priv_key
          }
        }

      private:
        const SchemeParameters params_;
        std::vector<AutomorphismDescription> alphas_;
        std::vector<AutomorphismDescription> betas_;
        std::vector<int> u_;
        std::vector<int> v_;
        std::vector<AutmorphismDescription> s_;
        std::vector<AutmorphismDescription> r_;

        AutomorphismDescription priv_key_base_;

        AutomorphismDescription c_;

        AutomorphismInvPair priv_key_;

        PublicKeys pub_keys_;

        AutomorphismInvPair shared_key_;

        static AutomorphismDescription get_automorphism_composition(const AutomorphismDescription& a1,
                                                                    const AutomorphismDescription& a2,
                                                                    const std::vector<int>& pattern) {
          auto pick = [&a1, &a2] (int i) {
            switch(i) {
              case -2:
                return a2.description_of_inverse();
              case -1:
                return a1.description_of_inverse();
              case 1:
                return a1;
              case 2:
                return a2;
            }
            assert(false);//should not get here
          };
          AutomorphismDescription result;
          for (auto i: pattern) {
            result *= pick(i);
          }
          return result;
        }

        //! Using public key calculate the private key.
        AutomorphismInvPair calculate_private_key(const PublicKeys& keys) {

        }

     };
  }
}

#endif // CRAG_FREE_GROUPS_FGACRYPTO_H
