/*
 * AAGCrypto.h
 *
 *  Created on: July 10, 2013
 *      Author: pmorar
 */

#ifndef CRAG_FREE_GROUPS_AAGCRYPTO_H
#define CRAG_FREE_GROUPS_AAGCRYPTO_H


#include <chrono>
#include "EndomorphismSLP.h"

namespace crag {

namespace aag_crypto {

  struct SchemeParameters {
      SchemeParameters(unsigned int rank,
                       unsigned int alice_tuple_size, unsigned int bob_tuple_size,
                       unsigned int lower_bound_pub_key_length, unsigned int upper_bound_pub_key_length,
                       unsigned int key_length)
        : RANK(rank), ALICE_TUPLE_SIZE(alice_tuple_size), BOB_TUPLE_SIZE(bob_tuple_size),
          LOWER_BOUND_PUB_KEY_LENGTH(lower_bound_pub_key_length), UPPER_BOUND_PUB_KEY_LENGTH(upper_bound_pub_key_length),
          KEY_LENGTH(key_length) {
        assert (LOWER_BOUND_PUB_KEY_LENGTH < UPPER_BOUND_PUB_KEY_LENGTH);
      }

      const unsigned int RANK;
      const unsigned int ALICE_TUPLE_SIZE;
      const unsigned int BOB_TUPLE_SIZE;
      const unsigned int LOWER_BOUND_PUB_KEY_LENGTH;
      const unsigned int UPPER_BOUND_PUB_KEY_LENGTH;
      const unsigned int KEY_LENGTH;
  };

  typedef EndomorphismSLP<int> Aut;
  typedef AutomorphismDescription<Aut > AutDescription;

  class PublicKey {
    public:
      PublicKey() = delete;

      //! Create random keys according to the specified parameters and random generator.
      template<typename RandomEngine, typename AutRandomGenerator>
      PublicKey(unsigned int tuple_size, RandomEngine& aut_size_generator, AutRandomGenerator& aut_rnd) {
        tuple_.reserve(tuple_size);
        for (unsigned int i = 0; i < tuple_size; ++i) {
          auto n = aut_size_generator();
          tuple_.push_back(AutDescription(n, aut_rnd));
        }
      }

      const AutDescription& operator[](std::size_t n) const {
        assert(n < num());
        return tuple_[n];
      }

      std::size_t num() const {
        return tuple_.size();
      }

    private:
      std::vector<AutDescription> tuple_;
  };

  class PrivateKey {
    public:
      PrivateKey() = delete;

      //! Create random keys according to the specified parameters and random generator.
      template<typename RandomEngine>
      PrivateKey(const PublicKey& pub_key, unsigned int length, RandomEngine& rand)
        : pub_key_(pub_key) {
        std::uniform_int_distribution<unsigned int> binary_dist(0, 1);
        std::uniform_int_distribution<std::size_t> pick_dist(0, pub_key.num() - 1);

        indices_.reserve(length);
        inverses_.reserve(length);

        for (unsigned int i = 0; i < length; ++i) {
          bool inverse = binary_dist(rand) == 1 ? true : false;
          inverses_.push_back(inverse);
          std::size_t n = pick_dist(rand);
          indices_.push_back(n);
          const auto& aut = pub_key[n];
          k_ *= inverse ? aut.inverse_description() : aut;
        }
      }

      const AutDescription& operator()() const {
        return k_;
      }

    private:
      const PublicKey& pub_key_;
      AutDescription k_;
      std::vector<std::size_t> indices_;
      std::vector<bool> inverses_;
      friend class KeysGenerator;
  };

  //! Provides blocks for ABA^{-1} and for BAB^{-1}
  class TransmittedInfo {
    public:
      TransmittedInfo() = delete;

      TransmittedInfo(const PublicKey& pub_key, const PrivateKey& pr_key) {
        tuple_.reserve(pub_key.num());
        for (std::size_t i = 0; i < pub_key.num(); ++i) {
          auto element = pub_key[i].conjugate_with(pr_key());
          tuple_.push_back(element.free_reduction().normal_form());
        }
      }

      const AutDescription& operator[](std::size_t n) const {
        assert(n < num());
        return tuple_[n];
      }

      std::size_t num() const {
        return tuple_.size();
      }

    private:
      std::vector<AutDescription> tuple_;
  };

  enum Mode {
    Alice,
    Bob
  };

  enum CalculationType {
    BlockFrNf,
    IterativeFrNf,
    SinlgeFrNf
  };

  typedef std::default_random_engine RandomEngine;

  //! Generates keys for the scheme.
  class KeysGenerator {
    public:

      KeysGenerator() = delete;

      //! Create random keys according to the specified parameters and random generator.
      KeysGenerator(const SchemeParameters& params, RandomEngine* p_rand)
        : params_(params),
          p_rand_(p_rand),
          length_distr_(params.LOWER_BOUND_PUB_KEY_LENGTH, params.UPPER_BOUND_PUB_KEY_LENGTH),
          aut_generator_(params.RANK, p_rand) {}

      PublicKey generate_public_key(Mode mode) {
        int tuple_size = mode == Alice ? params_.ALICE_TUPLE_SIZE : params_.BOB_TUPLE_SIZE;
        auto gen = std::bind(length_distr_, *p_rand_);
        return PublicKey(tuple_size, gen, aut_generator_);
      }

      PrivateKey generate_private_key(const PublicKey& pub_key) {
        return PrivateKey(pub_key, params_.KEY_LENGTH, *p_rand_);
      }



      Aut make_shared_key(const PrivateKey& priv_key, const TransmittedInfo& info, Mode mode, CalculationType calc_type) {
        //the shared key is aba^{-1}b^{-1}
        auto time = [] () {
            return std::chrono::high_resolution_clock::now();
          };
        auto start_time = time();
        auto duration = [&start_time, &time] () {
          return std::chrono::duration_cast<std::chrono::milliseconds>(time() - start_time);
        };

        if (mode == Alice) {
          Aut conj;
          auto inv_iter = priv_key.inverses_.crbegin();
          auto ind_iter = priv_key.indices_.crbegin();
          if (calc_type == BlockFrNf) {
            std::cout << "building blocks" << std::endl;
            unsigned int block_lengths = 5;
            std::vector<Aut> blocks;
            blocks.reserve(params_.KEY_LENGTH / block_lengths + 1);
            for (std::size_t i = 0; i < params_.KEY_LENGTH; i += block_lengths) {
              Aut block;
              for (std::size_t j = 0; j < block_lengths && j < params_.KEY_LENGTH - i; ++j) {
                auto k = i + j;
                bool inverse = priv_key.inverses_[k];
                auto n = priv_key.indices_[k];
                const auto& aut = info[n];
                block *= inverse ? aut() : aut.inverse();
              }

              std::cout << "block size " << slp_vertices_num(block);

              start_time = time();
              auto fr_block = block.free_reduction();
              std::cout << " fr_block size " << slp_vertices_num(fr_block)
                << " fr_block duration " << duration().count() << "ms";

              start_time = time();
              auto nf_block = fr_block.normal_form();
              std::cout << " nf_block size " << slp_vertices_num(nf_block)
                << " nf_block duration " << duration().count() << "ms" << std::endl;

              blocks.push_back(nf_block);
            }

            std::cout << "building key" << std::endl;
            for (const auto& block: blocks) {
              auto prod = conj * block;
              std::cout << " size " << slp_vertices_num(prod);
              start_time = time();
              auto fr = prod.free_reduction();
              std::cout << " fr size " << slp_vertices_num(fr)
                << " fr duration " << duration().count() << "ms";
              start_time = time();
              auto nf = fr.normal_form();
              std::cout << " nf size " << slp_vertices_num(nf)
                << " nf duration " << duration().count() << "ms" << std::endl;
              conj = nf;
            }
          } else if (calc_type == IterativeFrNf) {
            std::cout << "building key" << std::endl;
            for (; inv_iter != priv_key.inverses_.crend(); ++inv_iter, ++ind_iter) {
              bool inverse = *inv_iter;
              auto n = *ind_iter;
              const auto& aut = info[n];
              std::cout << "aut size " << slp_vertices_num(aut());
              auto prod = conj * (inverse ? aut() : aut.inverse());
              start_time = time();
              auto fr = prod.free_reduction();
              std::cout << " fr size " << slp_vertices_num(fr)
                << " fr duration " << duration().count() << "ms";
              start_time = time();
              auto nf = fr.normal_form();
              std::cout << " nf size " << slp_vertices_num(nf)
                << " nf duration " << duration().count() << "ms" << std::endl;
              conj = nf;
            }  
          } else {//SingleFrNf
            std::cout << "building key" << std::endl;
            for (; inv_iter != priv_key.inverses_.crend(); ++inv_iter, ++ind_iter) {
              bool inverse = *inv_iter;
              auto n = *ind_iter;
              const auto& aut = info[n];
              conj *= inverse ? aut() : aut.inverse();
            }
          }     
          
          auto key = priv_key()() * conj;

          std::cout << "size " << slp_vertices_num(key);

          start_time = time();  

          auto fr_key = key.free_reduction();
          std::cout << " fr size " << slp_vertices_num(fr_key)
            << " fr duration " << duration().count() << "ms";

          start_time = time();
          auto nf_key = fr_key.normal_form();
          std::cout << " nf size " << slp_vertices_num(nf_key)
            << " nf duration " << duration().count() << "ms" << std::endl;
          return nf_key;
        } else {//Bob
          Aut conj;
          auto inv_iter = priv_key.inverses_.cbegin();
          auto ind_iter = priv_key.indices_.cbegin();
          for (; inv_iter != priv_key.inverses_.cend(); ++inv_iter, ++ind_iter) {
            bool inverse = *inv_iter;
            auto n = *ind_iter;
            const auto& aut = info[n];
            conj *= inverse ? aut.inverse() : aut();
          }
          return conj * priv_key().inverse();
        }
      }



    private:
      const SchemeParameters params_;
      RandomEngine* p_rand_;
      std::uniform_int_distribution<typename RandomEngine::result_type> length_distr_;
      UniformAutomorphismSLPGenerator<int, RandomEngine> aut_generator_;

  };

  template<typename RandomEngine>
  KeysGenerator make_keys_generator(const SchemeParameters& params, RandomEngine* p_rand) {
    return KeysGenerator(params, p_rand);
  }

} // namespace aag_crypto
} // namespace crag

#endif // CRAG_FREE_GROUPS_AAGCRYPTO_H
