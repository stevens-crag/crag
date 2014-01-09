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


  typedef EndomorphismSLP Aut;
  typedef AutomorphismDescription<Aut> AutDescription;


  class PublicKey {
    public:
      PublicKey() = delete;

      //! Create random keys according to the specified parameters and random generator.
      template<typename RandomEngine, typename AutRandomGenerator>
      PublicKey(unsigned int key_num, RandomEngine& aut_size_generator, AutRandomGenerator& aut_rnd) {
        subkeys_.reserve(key_num);
        for (unsigned int i = 0; i < key_num; ++i) {
          subkeys_.push_back(AutDescription(aut_size_generator(), aut_rnd));
        }
      }

      const AutDescription& operator[](std::size_t n) const {
        assert(n < size());
        return subkeys_[n];
      }

      std::size_t size() const {
        return subkeys_.size();
      }

    private:
      std::vector<AutDescription> subkeys_;
  };


  class PrivateKey {
    public:

      struct PartDescription {
        const std::size_t public_key_index;
        const bool inversed;

        PartDescription(std::size_t public_key_index, bool inversed)
          : public_key_index(public_key_index),
            inversed(inversed) {}
      };

      PrivateKey() = delete;

      //! Create random keys according to the specified parameters and random generator.
      template<typename RandomEngine>
      PrivateKey(const PublicKey& pub_key, std::size_t length, RandomEngine& rand)
        : pub_key_(pub_key),
          k_() {
        std::uniform_int_distribution<std::size_t> binary_distribution(0, 1);
        std::uniform_int_distribution<std::size_t> discrete_distribution(0, pub_key.size() - 1);

        indices_.reserve(length);
        inverses_.reserve(length);

        for (std::size_t i = 0; i < length; ++i) {
          bool is_inverse = binary_distribution(rand) == 1;
          inverses_.push_back(is_inverse);
          std::size_t n = discrete_distribution(rand);
          indices_.push_back(n);
          const auto& aut = pub_key[n];
          k_ *= is_inverse ? aut.inverse_description() : aut;
        }
      }

      const AutDescription& operator()() const {
        return k_;
      }

      //! Returns the length in terms of public keys composing it.
      std::size_t length() const {
        return indices_.size();
      }

      //! Returns the description of the public key which is the nth part of the private key.
      PartDescription part_description(std::size_t n) const {
        return PartDescription(indices_[n], inverses_[n]);
      }

    private:
      const PublicKey& pub_key_;
      AutDescription k_;
      std::vector<std::size_t> indices_;
      std::vector<bool> inverses_;
  };


  //! Provides building blocks for ABA^{-1} and for BAB^{-1}, which parties send to each other.
  class TransmittedInfo {
    public:
      TransmittedInfo() = delete;

      TransmittedInfo(const PublicKey& public_key, const PrivateKey& private_key) {
        subkeys_.reserve(public_key.size());
        for (std::size_t i = 0; i < public_key.size(); ++i) {
          auto element = public_key[i].conjugate_with(private_key());
          subkeys_.push_back(element.free_reduction().normal_form());
        }
      }

      const AutDescription& operator[](std::size_t n) const {
        assert(n < size());
        return subkeys_[n];
      }

      std::size_t size() const {
        return subkeys_.size();
      }

    private:
      std::vector<AutDescription> subkeys_;
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
  template<typename RandomEngine = std::default_random_engine>
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
        auto get_part = [&priv_key, &info, mode] (std::size_t index) {
          if (mode == Participant::Alice) {
            auto part = priv_key.part_description(priv_key.length() - 1 - index);
            const auto& aut = info[part.public_key_index];
            return part.inversed ? aut() : aut.inverse();
          } else {
            auto part = priv_key.part_description(index);
            const auto& aut = info[part.public_key_index];
            return part.inversed ? aut.inverse() : aut();
          }
        };

        if (mode == Participant::Alice) {
          std::cout << "Alice is building key" << std::endl;
        } else {
          std::cout << "Bob is building key" << std::endl;
        }

        Aut conjugation;
        switch (calc_type) {
          case CalculationType::BlockReduction: {
              std::cout << "BlockReduction" << std::endl;
              std::vector<Aut> blocks;
              blocks.reserve(params_.KEY_LENGTH / block_length + 1);
              std::size_t index = 0;
              for (std::size_t i = 0; i < params_.KEY_LENGTH; i += block_length) {
                Aut block;
                for (std::size_t j = 0; j < block_length && j < params_.KEY_LENGTH - i; ++j) {
                  block *= get_part(index++);
                }
                blocks.push_back(AutomorphismReducer::reduce(block, true));
              }

              std::cout << "building key" << std::endl;
              for (const auto& block: blocks) {
                conjugation = AutomorphismReducer::reduce(conjugation * block, true);
              }
              break;
            }
          case CalculationType::IterativeReduction: {
              std::cout << "IterativeReduction" << std::endl;
              for (std::size_t index = 0; index < priv_key.length(); ++index) {
                conjugation = AutomorphismReducer::reduce(conjugation * get_part(index), true);
              }
              break;
            }
          case CalculationType::ThresholdReduction: {
              std::cout << "ThresholdReduction" << std::endl;
              for (std::size_t index = 0; index < priv_key.length(); ++index) {
                conjugation *= get_part(index);
                auto vertices_num = slp_vertices_num(conjugation);
                std::cout << "key part |" << vertices_num << "|" << std::endl;
                if (vertices_num > fold_threshold) {
                  std::cout << "size exceeded threshold: folding..." << std::endl;
                  conjugation = AutomorphismReducer::reduce(conjugation, true);
                }
              }
              break;
            }
          case CalculationType::SingleReduction: {
              std::cout << "SingleReduction" << std::endl;
              for (std::size_t index = 0; index < priv_key.length(); ++index) {
                conjugation *= get_part(index);
              }
              break;
            }
          default:
            throw std::invalid_argument("Invalid argument!");
        }

        std::cout << "key ";
        if (mode == Participant::Alice) {
          return AutomorphismReducer::reduce(priv_key().aut() * conjugation, true);
        } else {//Bob
          return AutomorphismReducer::reduce(conjugation * priv_key().inverse(), true);
        }
      }

    private:
      const SchemeParameters params_;
      RandomEngine* p_rand_;
      std::uniform_int_distribution<typename RandomEngine::result_type> length_distr_;
      UniformAutomorphismSLPGenerator<RandomEngine> aut_generator_;

//      template<typename Iterator>
//      Aut combine_keys()

  };

  template<typename RandomEngine>
  KeysGenerator<RandomEngine> make_keys_generator(const SchemeParameters& params, RandomEngine* p_rand) {
    return KeysGenerator<RandomEngine>(params, p_rand);
  }

} // namespace aag_crypto
} // namespace crag

#endif // CRAG_FREE_GROUPS_AAGCRYPTO_H
