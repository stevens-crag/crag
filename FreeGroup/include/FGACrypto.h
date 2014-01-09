/*
 * FGACrypto.h
 *
 *  Created on: Apr 27, 2013
 *      Author: pmorar
 */

#ifndef CRAG_FREE_GROUPS_FGACRYPTO_H
#define CRAG_FREE_GROUPS_FGACRYPTO_H

#include <chrono>
#include "EndomorphismSLP.h"

#include <fstream>

/**
 * The file contains cryptoscheme programmed corresponding to the specifications by Sasha Ushakov.
 */
namespace crag {

namespace fga_crypto {

// uncomment for debugging info
//#define CRAG_FGA_CRYPTO_DEBUG_OUTPUT

  struct SchemeParameters {

      SchemeParameters(int n, int u, int v, int c)
        : N(n), A_SIZE(5), B_SIZE(4), U_LENGTH(u), V_LENGTH(v), C_SIZE(c) {}

      SchemeParameters(int n, int a, int b, int u, int v, int c)
        : N(n), A_SIZE(a), B_SIZE(b), U_LENGTH(u), V_LENGTH(v), C_SIZE(c) {}


      const unsigned int N = 3;
      const unsigned int A_SIZE = 5;
      const unsigned int B_SIZE = 4;
      const unsigned int U_LENGTH = 25;
      const unsigned int V_LENGTH = 25;
      const unsigned int C_SIZE = 50;
  };

  typedef EndomorphismSLP Aut;
  typedef AutomorphismDescription<Aut> AutDescription;
  typedef CommutatorSet<AutDescription> CommSet;


  //! Public keys
  class PublicKeys {
    public:

      PublicKeys() {}

      PublicKeys(const std::vector<AutDescription>& alphas, const std::vector<AutDescription>& betas) {
        s_.reserve(4);
        for (int i = 0; i < 4; ++i) {
          s_.push_back(alphas[i] * betas[i]);
#ifdef CRAG_FGA_CRYPTO_DEBUG_OUTPUT
          std::cout << "s " << i + 1 << " n=" << s_[i].composed_num() <<  std::endl;
          s_[i]().print(&std::cout);
#endif
        }

        r_.reserve(4);
        for (int i: {0, 1}) {
          for (int j: {2, 3}) {
            r_.push_back(CommSet(betas[i], betas[j]));
#ifdef CRAG_FGA_CRYPTO_DEBUG_OUTPUT
            std::cout << "r_" << i + 1 << "_" << j + 1 << " n=" << r_.back().get(false, false).composed_num() <<  std::endl;
            r_.back().get(false, false)().print(&std::cout);
#endif
          }
        }
      }

      AutDescription get_s(int index) const {
        if (index > 0) {
          return s_[index - 1];
        } else {
          return s_[- index - 1].inverse_description();
        }
      }

      const AutDescription& get_r(int first_index, int second_index) const {
        bool first_inversed = first_index < 0;
        bool second_inversed = second_index < 0;
        first_index = first_inversed ? - first_index : first_index;
        first_index -= 1;
        second_index = second_inversed ? - second_index : second_index;
        second_index -= 3;
        return r_[2 * first_index + second_index].get(first_inversed, second_inversed);
      }

      //! Make the normal form of free reduction
      PublicKeys reduce() const {
        std::vector<AutDescription> new_s;
        new_s.reserve(s_.size());
        std::transform(s_.cbegin(), s_.cend(),
                       std::back_inserter(new_s),
                       [] (const AutDescription& ad) {return AutomorphismReducer::reduce(ad);});

        std::vector<CommSet> new_r;
        new_r.reserve(r_.size());
        std::transform(r_.cbegin(), r_.cend(),
                       std::back_inserter(new_r),
                       [] (const CommSet& cs) {return cs.reduce();});
        return PublicKeys(std::move(new_s), std::move(new_r));
      }

      PublicKeys conjugate_with(const AutDescription& conjugator) const {
        return PublicKeys(std::move(conjugate_all(s_, conjugator)),
                          std::move(conjugate_all(r_, conjugator)));
      }

    private:
      PublicKeys(std::vector<AutDescription>&& v_s, std::vector<CommSet>&& v_r)
        : s_(v_s), r_(v_r) {}

      std::vector<AutDescription> s_;
      std::vector<CommSet> r_;
  };

  enum class Mode {
      Alice,
      Bob
  };

  //! Generates keys for the scheme.
  class KeysGenerator {
    public:

      bool is_logging = false;

      KeysGenerator() = delete;

      //! Create random keys according to the specified parameters and random generator.
      template<typename RandomEngine>
      KeysGenerator(const SchemeParameters& params, RandomEngine& rand)
        : params_(params),
          alphas_(),
          betas_(),
          u_(),
          v_(),
          priv_key_base_(),
          c_(),
          priv_key_(),
          pub_keys_(),
          shared_key_() {

        //generating alphas, betas
        UniformAutomorphismSLPGenerator<RandomEngine> random_for_alphas(params.N, &rand);
        UniformAutomorphismSLPGenerator<RandomEngine> random_for_betas(params.N + 1, 2 * params.N, &rand);
        UniformAutomorphismSLPGenerator<RandomEngine> random_for_c(2 * params.N, &rand);

        alphas_.reserve(4);
        betas_.reserve(4);
        for (int i = 0; i < 4; ++i) {
          alphas_.push_back(AutDescription(params.A_SIZE, random_for_alphas));
          betas_.push_back(AutDescription(params.B_SIZE, random_for_betas));
        }

#ifdef CRAG_FGA_CRYPTO_DEBUG_OUTPUT
        for (int i = 0; i < 4; ++i) {
          std::cout << "alpha " << i + 1 << " n=" << alphas_[i].composed_num() << std::endl;
          alphas_[i]().print(&std::cout);
        }

        for (int i = 0; i < 4; ++i) {
          std::cout << "beta " << i + 1 << " n=" << betas_[i].composed_num() <<  std::endl;
          betas_[i]().print(&std::cout);
        }
#endif

        c_ = AutDescription(params.C_SIZE, random_for_c);

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
#ifdef CRAG_FGA_CRYPTO_DEBUG_OUTPUT
        std::cout << "u ";
        for (int i: u_)
          std::cout << i << " ";
        std::cout << std::endl;

        std::cout << "v ";
        for (int i: v_)
          std::cout << i << " ";
        std::cout << std::endl;
#endif
        pub_keys_ = PublicKeys(alphas_, betas_).conjugate_with(c_).reduce();
#ifdef CRAG_FGA_CRYPTO_DEBUG_OUTPUT
        std::cout << "pub key" << std::endl;
        std::cout << "s" << std::endl;
        for (auto& ad: pub_keys_.s) {
          AutomorphismReducer::reduce(ad, true);
        }
        std::cout << "r" << std::endl;
        for (auto& cs: pub_keys_.r) {
          AutomorphismReducer::reduce(cs.get(false, false), true);
        }
#endif

        priv_key_base_ = make_commutator(get_betas_composition(u_),
            get_betas_composition(v_));
        priv_key_ = priv_key_base_.conjugate_with(c_);
#ifdef CRAG_FGA_CRYPTO_DEBUG_OUTPUT
        std::cout << "priv key base" << " n=" << priv_key_base_.composed_num() <<   std::endl;
        AutomorphismReducer::reduce(priv_key_base_, true);
        std::cout << "priv key" << " n=" << priv_key_.composed_num() <<  std::endl;
        AutomorphismReducer::reduce(priv_key_, true);
#endif
      }

      const PublicKeys& public_keys() const {
        return pub_keys_;
      }

      const AutDescription& private_key() const {
        return priv_key_;
      }

      //! Processes public keys provided by other party to send them back
      PublicKeys process_incoming_public_keys(const PublicKeys& incoming_public_keys) {
        return incoming_public_keys.conjugate_with(priv_key_).reduce();
      }

      //! Make shared keys with the given public keys of another party
      /**
       * @param processed_public_keys our public keys processed by another party
       * @param order if true Alice makes the key, otherwise Bob
       * @return
       */
      Aut make_shared_key(const PublicKeys& processed_public_keys, Mode mode) {
        Aut key;
        std::map<int, Aut> line_cache;
        const std::vector<int>& indices = mode == Mode::Alice ? u_ : v_;
        for (int i = indices.size() - 1; i >= 0; --i) {
          const int row_index = indices[i];
          auto cached_ad = line_cache.find(row_index);
          Aut line;
          if (cached_ad != line_cache.end()) {
            line = cached_ad->second;
          } else {
            if (is_logging) {
              std::cout << "new line calculation:" << std::endl;
            }
            line = calculate_private_key_line(row_index, processed_public_keys, mode);
            line_cache.insert(std::make_pair(row_index, line));
          }
          if (is_logging) {
            std::cout << "mult:" << std::endl;
          }
          key = line * key;
          if (i > 0) {
            const int conj_index = indices[i - 1];
            key = key.conjugate_with(processed_public_keys.get_s(conj_index));
          }
          if (is_logging) {
            std::cout << "line addtion: ";
          }
          key = AutomorphismReducer::reduce(key, is_logging);
        }

        if (mode == Mode::Alice) {
          key = priv_key_() * key;//Alice: a * (bab^-1)^-1
        } else {
          key *= priv_key_.inverse();//Bob: aba^-1 *= b^-1
        }

        if (is_logging) {
          std::cout << "key: ";
        }
        key = AutomorphismReducer::reduce(key, is_logging);
        return key;
      }

    private:
      const SchemeParameters params_;
      std::vector<AutDescription> alphas_;
      std::vector<AutDescription> betas_;
      std::vector<int> u_;
      std::vector<int> v_;

      AutDescription priv_key_base_;

      AutDescription c_;

      AutDescription priv_key_;

      PublicKeys pub_keys_;

      AutDescription shared_key_;

      AutDescription get_betas_composition(const std::vector<int>& pattern) {
        auto pick = [&] (int i) {
          if (i > 0) {
            return betas_[i - 1];
          } else {
            return betas_[- i - 1].inverse_description();
          }
        };
        AutDescription result;
        for (auto i: pattern) {
          result *= pick(i);
        }
        return result;
      }

      Aut calculate_private_key_line(int row_index, const PublicKeys& public_key, Mode mode) {
        Aut value;
        const std::vector<int>& indices = mode == Mode::Alice ? v_ : u_;
        auto get_r = [&mode, &row_index, &public_key] (int col_index) {
          return mode == Mode::Alice
              ? public_key.get_r(row_index, col_index).inverse()
              : public_key.get_r(col_index, row_index)();
        };

        for (int i = indices.size() - 1; i >= 0; --i) {
          const int col_index = indices[i];
          value *= get_r(col_index);
          const int conj_index = indices[i-1];
          if (i > 0) {
            value = value.conjugate_with(public_key.get_s(conj_index));
          }
          if (is_logging) {
            std::cout << "line step ";
          }
          value = AutomorphismReducer::reduce(value, is_logging);//reducing the size
        }
        return value;
      }
   };
} // namespace fga_crypto
} // namespace crag

#endif // CRAG_FREE_GROUPS_FGACRYPTO_H
