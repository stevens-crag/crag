#include "walnut_attack.h"

using namespace crag;

using ZZ5 = finitefield::ZZ<5>;
using GF32 = finitefield::FieldElement<finitefield::IdealGeneratedByPolynomial<finitefield::ZZ<2>, 1, 0, 1, 0, 0, 1>>;
using GF256 =
    finitefield::FieldElement<finitefield::IdealGeneratedByPolynomial<finitefield::ZZ<2>, 1, 1, 0, 1, 1, 0, 0, 0, 1>>;


static void printTotalxNumbersSeq3(int N, const Word& w, size_t strand_a, size_t strand_b) {
  {
    int total = 0;
    std::vector<std::vector<int>> result(N, vector<int>(N, 0));
    Permutation p1(N);
    Permutation p2(N);


    int i = 0;
    std::cout << "Show[" << std::endl;
    std::cout << "ListPlot[{";
    for (const auto let : w) {
      const auto ind = std::abs(let);

      auto a = p2[ind - 1];
      auto b = p2[ind];
      p1.change(a, b);
      p2.change(ind - 1, ind);
      if (a < b) {
        std::swap(a, b);
      }

      if (i != 0) {
        std::cout << ",";
      }

      std::cout << "{" << i << "," << total << "}";

      const auto old_val = std::abs(result[a][b]);

      result[a][b] += (let < 0 ? -1 : 1);

      if (old_val < std::abs(result[a][b])) {
        total++;
      } else {
        total--;
      }
      ++i;
    }
    std::cout << "}, PlotStyle -> LightGray], " << std::endl;
  }

  {
    int total = 0;
    std::vector<std::vector<int>> result(N, vector<int>(N, 0));
    Permutation p1(N);
    Permutation p2(N);

    int i = 0;
    std::cout << "ListPlot[{";
    for (const auto let : w) {
      const auto ind = std::abs(let);

      auto a = p2[ind - 1];
      auto b = p2[ind];
      p1.change(a, b);
      p2.change(ind - 1, ind);
      if (a < b) {
        std::swap(a, b);
      }

      if (strand_a == b && strand_b == a) {
        std::cout << ",";
        std::cout << "{" << i << "," << total << "}";
      }

      const auto old_val = std::abs(result[a][b]);

      result[a][b] += (let < 0 ? -1 : 1);

      if (old_val < std::abs(result[a][b])) {
        total++;
      } else {
        total--;
      }
      ++i;
    }

    std::cout << "}, PlotStyle -> Red]" << std::endl;
    std::cout << "]" << std::endl;
  }
}

int main() {
  size_t success_count = 0;

  const size_t init_seed = 0;
  const size_t experiments_count = 100;

  for (size_t e = init_seed; e < init_seed + experiments_count; ++e) {
    std::cout << "========================================================" << std::endl;
    std::cout << "Experiment #" << e << std::endl;

    // get random protocol instance with seed = e
    const auto protocol = walnut::getProtocolFor128BitsSecurity<walnut::StabilizerSquare>(e);
    // const auto protocol = walnut::getProtocolFor128BitsSecurity<walnut::StabilizerDoubleSquare>(e);

    // const auto protocol = walnut::getProtocolFor256BitsSecurity<walnut::StabilizerSquare>(e);
    // const auto protocol = walnut::getProtocolFor256BitsSecurity<walnut::StabilizerDoubleSquare>(e);

    //    const auto protocol = walnut::getProtocolFor256BitsSecurityN11<walnut::StabilizerSquare>(e);
    //    const auto protocol = walnut::getProtocolFor256BitsSecurityN11<walnut::StabilizerDoubleSquare>(e);

    const auto n = protocol.publicParameters().n();

    // generate random private key using seed = e
    const auto private_key = protocol.generatePrivateKey(e);
    const auto public_key = protocol.computePublicKey(private_key);

    std::cout << "|w1|  = " << private_key.w1().length() << ", ";
    std::cout << "|w2|  = " << private_key.w2().length() << std::endl;

    // We have several attempts to find the keys
    int spent_time = time(0);
    bool success = false;
    const size_t attempt_number = 10;

    for (size_t attempt = 0u; (attempt < attempt_number) && !success; ++attempt) {
      std::cout << "Attempt #" << attempt << std::endl;

      std::mt19937_64 g(attempt);

      const auto fake_private_key = walnut::attack(protocol, private_key, public_key, g);

      if (success = (fake_private_key != boost::none)) {
        std::cout << "|fake w1|  = " << shortenBraid2(n, fake_private_key->w1()).length() << ", ";
        std::cout << "|fake w2|  = " << shortenBraid2(n, fake_private_key->w2()).length() << std::endl;

        if (success &= walnut::checkFakePrivateKey(protocol, *fake_private_key, private_key, public_key, g)) {
          std::cout << "Computed private key is GOOD." << std::endl;
        } else {
          std::cout << "Computed private key is BAD." << std::endl;
        }
      }
    }

    spent_time = time(0) - spent_time;

    if (success) {
      success_count++;
      std::cout << "Experiment #" << e << " is successful" << std::endl;
    } else {
      std::cout << "Experiment #" << e << " failed" << std::endl;
    }

    std::cout << "   Time = " << spent_time << std::endl;
    std::cout << "Success: " << success_count << " out of " << (e + 1 - init_seed) << std::endl;
  }

  return 0;
}
