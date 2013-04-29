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
        const unsigned int C_LENGTH = 25;
        //TODO make constructor
    };

    /**
     * Automorphism description containing parts of which it is composed and its inverse.
     * The class is effectively immutable.
     */
    template<typename TerminalSymbol = int>
    class AutomorphismDescription {
      public:
        template<typename RandomAutomorphismGenerator>
        AutomorphismDescription(unsigned int size, RandomAutomorphismGenerator& random) {
          a_parts_.reserve(size);
          a_inv_parts_.reserve(size);
          for (unsigned int i = 0; i < size; ++i) {
            a_parts_.push_back(random());
          }
          std::for_each(a_parts_.rbegin(), a_parts_.rend(), [&a_inv_parts_] (const EndomorphismSLP<TerminalSymbol>& e) {
            a_inv_parts_.push_back(e.inverse());
          });
          a_ = EndomorphismSLP<TerminalSymbol>::composition(a_parts_.begin(), a_parts_.end());
          a_inv_ = EndomorphismSLP<TerminalSymbol>::composition(a_inv_parts_.begin(), a_inv_parts_.end());
        }

        //! Returns the automorphism itself.
        const EndomorphismSLP<TerminalSymbol>& operator()() const {
          return a_;
        }

        //! Returns the automorphism inverse.
        const EndomorphismSLP<TerminalSymbol>& inverse() const {
          return a_inv_;
        }

        //! Apply the provided function to each of the composition parts of the automorphism.
        template<typename Func>
        Func for_each_comp_part(Func f) const {
          return std::for_each(a_parts_.begin(), a_parts_.end(), f);
        }

        //! Apply the provided function to each of the composition parts of the autmorphism inverse.
        template<typename Func>
        Func for_each_comp_part(Func f) const {
          return std::for_each(a_inv_parts_.begin(), a_inv_parts_.end(), f);
        }

      private:
        EndomorphismSLP<TerminalSymbol> a_;
        EndomorphismSLP<TerminalSymbol> a_inv_;
        std::vector<EndomorphismSLP<TerminalSymbol> > a_parts_;
        std::vector<EndomorphismSLP<TerminalSymbol> > a_inv_parts_;
    };

    class Keys {
      public:

        Keys() = delete;

        template<typename RandomEngine>
        Keys(const SchemeParameters& params, RandomEngine& rand)
          : params_(params),
            alphas_(),
            betas_(),
            u_(),
            v_(),
            c_() {
          //generating keys

          //generating alphas and betas
          UniformAutomorphismSLPGenerator<int, RandomEngine> random(params.N, p_rand);
          alphas_.reserve(4);
          betas_.reserve(4);
          for (int i = 0; i < 4; ++i) {
            alphas_.push_back(EndomorphismSLP<int>::composition(params.A_SIZE, random));
            betas_.push_back(EndomorphismSLP<int>::composition(params.B_SIZE, random));
          }

          //generating u, v

          //generating c

        }

      private:
        const SchemeParameters params_;
        std::vector<EndomorphismSLP<int> > alphas_;
        std::vector<EndomorphismSLP<int> > alphas_;
        std::vector<EndomorphismSLP<int> > betas_;
        std::vector<int> u_;
        std::vector<int> v_;
        EndomorphismSLP<int> c_;

     };
  }
}

#endif // CRAG_FREE_GROUPS_FGACRYPTO_H
