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
        const unsigned int C_SIZE = 59;
        //TODO make constructor
    };

    struct PublicKeys {
        std::vector<AutmorphismDescription> s;
        std::vector<AutmorphismDescription> r;
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
            c_(params.C_SIZE, UniformAutomorphismSLPGenerator<int, RandomEngine>(2 * params.N, p_rand)),
            priv_key_(),
            pub_keys_() {
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

          //generating public keys
          pub_keys_.s.reserve(4);
          for (int i = 0; i < 4; ++i) {
            pub_keys_.s.push_back(alphas_[i] * betas_[i]);
          }

          auto commutator = [] (const AutomorphismDescription& a1, const AutomorphismDescription& a2) {
            return a1 * a2 * a1.invert() * a2.invert();
          };

          pub_keys_.r.reserve(4);
          for (unsigned int i = 0; i <= 1; ++i)
            for (unsigned int j = 3; j <= 4; ++j)
              pub_keys_.r.push_back(commutator(betas_[i], betas_[j]));

          //generating private key
          priv_key_ = commutator(get_automorphism_composition(betas_[0], betas_[1], u_),
              get_automorphism_composition(betas_[2], betas_[3], v_));

        }

      private:
        const SchemeParameters params_;
        std::vector<AutomorphismDescription> alphas_;
        std::vector<AutomorphismDescription> betas_;
        std::vector<int> u_;
        std::vector<int> v_;
        AutomorphismDescription c_;

        AutmorphismDescription priv_key_;

        PublicKeys pub_keys_;

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

     };


    template
    AutomorphismDescription invert() const {

    }
  }
}

#endif // CRAG_FREE_GROUPS_FGACRYPTO_H
