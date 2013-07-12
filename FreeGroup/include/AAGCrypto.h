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

  const std::size_t fold_threshold = 50000;

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
  typedef std::default_random_engine RandomEngine;


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
      PrivateKey(const PublicKey& pub_key, std::size_t length, RandomEngine& rand)
        : pub_key_(pub_key),
          k_() {
        std::uniform_int_distribution<std::size_t> binary_dist(0, 1);
        std::uniform_int_distribution<std::size_t> pick_dist(0, pub_key.num() - 1);

        indices_.reserve(length);
        inverses_.reserve(length);

        for (std::size_t i = 0; i < length; ++i) {
          bool is_inverse = binary_dist(rand) == 1;
          inverses_.push_back(is_inverse);
          std::size_t n = pick_dist(rand);
          indices_.push_back(n);
          const auto& aut = pub_key[n];
          k_ *= is_inverse ? aut.inverse_description() : aut;
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


  //! Provides building blocks for ABA^{-1} and for BAB^{-1}, which parties send to each other.
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


  enum class Participant {
    Alice,
    Bob
  };


  enum class CalculationType {
    BlockReduction,
    IterativeReduction,
    ThresholdReduction,//use threshold
    SingleReduction
  };


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

      PublicKey generate_public_key(Participant mode) {
        int tuple_size = mode == Participant::Alice ? params_.ALICE_TUPLE_SIZE : params_.BOB_TUPLE_SIZE;
        auto gen = std::bind(length_distr_, *p_rand_);
        return PublicKey(tuple_size, gen, aut_generator_);
      }

      PrivateKey generate_private_key(const PublicKey& pub_key) {
        return PrivateKey(pub_key, params_.KEY_LENGTH, *p_rand_);
      }

      Aut make_shared_key(const PrivateKey& priv_key, const TransmittedInfo& info, Participant mode,
                          CalculationType calc_type, std::size_t block_length = 5) const {
        //the shared key is aba^{-1}b^{-1}
//        auto time = [] () {
//            return std::chrono::high_resolution_clock::now();
//          };
//        auto start_time = time();
//        auto duration = [&start_time, &time] () {
//          return std::chrono::duration_cast<std::chrono::milliseconds>(time() - start_time);
//        };
//        auto duration_in_ms = duration();

        if (mode == Participant::Alice) {
          Aut conj;
          auto inverses_iterator = priv_key.inverses_.crbegin();
          auto indices_iterator = priv_key.indices_.crbegin();
          switch (calc_type) {
            case CalculationType::BlockReduction: {
                std::cout << "building blocks" << std::endl;
                std::vector<Aut> blocks;
                blocks.reserve(params_.KEY_LENGTH / block_length + 1);
                for (std::size_t i = 0; i < params_.KEY_LENGTH; i += block_length) {
                  Aut block;
                  for (std::size_t j = 0;
                       j < block_length && j < params_.KEY_LENGTH - i;
                       ++j, ++inverses_iterator, ++indices_iterator) {
                    bool inverse = *inverses_iterator;
                    auto n = *indices_iterator;
                    const auto& aut = info[n];
                    block *= inverse ? aut() : aut.inverse();
                  }

                  std::cout << "block ";

                  blocks.push_back(AutomorphismReducer::reduce(block, true));
                }

                std::cout << "building key" << std::endl;
                for (const auto& block: blocks) {
                  auto prod = conj * block;
                  std::cout << "key part";
                  conj = AutomorphismReducer::reduce(prod, true);
                }
                break;
              }
            case CalculationType::IterativeReduction: {
                std::cout << "building key" << std::endl;
                for (; inverses_iterator != priv_key.inverses_.crend(); ++inverses_iterator, ++indices_iterator) {
                  bool inverse = *inverses_iterator;
                  auto n = *indices_iterator;
                  const auto& aut = info[n];
                  auto prod = conj * (inverse ? aut() : aut.inverse());
                  std::cout << "key part";
                  conj = AutomorphismReducer::reduce(prod, true);
                }
                break;
              }
            case CalculationType::ThresholdReduction: {
                std::cout << "building key" << std::endl;
                for (; inverses_iterator != priv_key.inverses_.crend(); ++inverses_iterator, ++indices_iterator) {
                  bool inverse = *inverses_iterator;
                  auto n = *indices_iterator;
                  const auto& aut = info[n];
                  conj *= inverse ? aut() : aut.inverse();
                  const std::size_t v_num = slp_vertices_num(conj);
                  std::cout << "key part |" << v_num << "|" << std::endl;
                  if (v_num > fold_threshold) {
                    std::cout << "size exceeded threshold: folding..." << std::endl;
                    conj = AutomorphismReducer::reduce(conj, true);
                  }
                }
                break;
              }
            case CalculationType::SingleReduction: {
                std::cout << "building key" << std::endl;
                for (; inverses_iterator != priv_key.inverses_.crend(); ++inverses_iterator, ++indices_iterator) {
                  bool inverse = *inverses_iterator;
                  auto n = *indices_iterator;
                  const auto& aut = info[n];
                  conj *= inverse ? aut() : aut.inverse();
                }
                break;
              }
            default:
              throw std::invalid_argument("Invalid argument!");
          }

          auto key = priv_key()() * conj;

          std::cout << "key ";
          return AutomorphismReducer::reduce(key, true);
        } else {//Bob
          assert (mode == Participant::Bob);
          Aut conj;
          auto inverses_iterator = priv_key.inverses_.cbegin();
          auto indices_iterator = priv_key.indices_.cbegin();
          switch (calc_type) {
            case CalculationType::BlockReduction: {
                std::cout << "building blocks" << std::endl;
                std::vector<Aut> blocks;
                blocks.reserve(params_.KEY_LENGTH / block_length + 1);
                for (std::size_t i = 0; i < params_.KEY_LENGTH; i += block_length) {
                  Aut block;
                  for (std::size_t j = 0;
                       j < block_length && j < params_.KEY_LENGTH - i;
                       ++j, ++inverses_iterator, ++indices_iterator) {
                    bool inverse = *inverses_iterator;
                    auto n = *indices_iterator;
                    const auto& aut = info[n];
                    block *= inverse ? aut.inverse() : aut();
                  }
                  std::cout << "block ";

                  blocks.push_back(AutomorphismReducer::reduce(block, true));
                }

                std::cout << "building key" << std::endl;
                for (const auto& block: blocks) {
                  auto prod = conj * block;
                  std::cout << "key part";
                  conj = AutomorphismReducer::reduce(prod, true);
                }
                break;
              }
            case CalculationType::IterativeReduction: {
                std::cout << "building key" << std::endl;
                for (; inverses_iterator != priv_key.inverses_.cend(); ++inverses_iterator, ++indices_iterator) {
                  bool inverse = *inverses_iterator;
                  auto n = *indices_iterator;
                  const auto& aut = info[n];
                  auto prod = conj * (inverse ? aut.inverse() : aut());
                  std::cout << "key part";
                  conj = AutomorphismReducer::reduce(prod, true);
                }
                break;
              }
            case CalculationType::ThresholdReduction: {
                std::cout << "building key" << std::endl;
                for (; inverses_iterator != priv_key.inverses_.cend(); ++inverses_iterator, ++indices_iterator) {
                  bool inverse = *inverses_iterator;
                  auto n = *indices_iterator;
                  const auto& aut = info[n];
                  conj *= inverse ? aut.inverse() : aut();
                  const std::size_t v_num = slp_vertices_num(conj);
                  std::cout << "key part |" << v_num << "|" << std::endl;
                  if (v_num > fold_threshold) {
                    std::cout << "size exceeded threshold: folding..." << std::endl;
                    conj = AutomorphismReducer::reduce(conj, true);
                  }
                }
                break;
              }
            case CalculationType::SingleReduction: {
                std::cout << "building key" << std::endl;
                for (; inverses_iterator != priv_key.inverses_.cend(); ++inverses_iterator, ++indices_iterator) {
                  bool inverse = *inverses_iterator;
                  auto n = *indices_iterator;
                  const auto& aut = info[n];
                  conj *= inverse ? aut.inverse() : aut();
                }
                break;
              }
            default:
              throw std::invalid_argument("Invalid argument!");
          }
          auto key = conj * priv_key().inverse();

          std::cout << "key ";
          return AutomorphismReducer::reduce(key, true);
        }
      }

    private:
      const SchemeParameters params_;
      RandomEngine* p_rand_;
      std::uniform_int_distribution<typename RandomEngine::result_type> length_distr_;
      UniformAutomorphismSLPGenerator<int, RandomEngine> aut_generator_;

//      template<typename Iterator>
//      Aut combine_keys()

  };


  template<typename RandomEngine>
  KeysGenerator make_keys_generator(const SchemeParameters& params, RandomEngine* p_rand) {
    return KeysGenerator(params, p_rand);
  }

} // namespace aag_crypto
} // namespace crag

#endif // CRAG_FREE_GROUPS_AAGCRYPTO_H
