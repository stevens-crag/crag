#include <fstream>

#include "ThLeftNormalForm.h"
#include "braid_group.h"
#include "kayawood_attack.h"

using namespace crag;

int main() {
  const size_t init_seed = 0;
  const size_t experiments_count = 100;

  size_t success_count = 0;

  for (size_t i = init_seed; i < init_seed + experiments_count; ++i) {
    std::cout << "Generating random instance with seed = " << i << std::endl;
//    const auto protocol = kayawood::getProtocolFor128BitsSecurity(i);
//    const auto protocol = kayawood::getProtocolFor256BitsSecurity(i);
//    const auto protocol = kayawood::getProtocolFor128BitsSecurityMultipleCloaking(i);
    const auto protocol = kayawood::getProtocolFor256BitsSecurityMultipleCloaking(i);

    const auto parameters = protocol.parameters();
    const auto instance = protocol.generateInstance(i);

    if (kayawood::isBadInstance(instance)) {
      std::cout << "Protocol instance #" << i << " is bad, shared key can be obtained from public keys." << std::endl;
    }

    const auto priv_a = instance.alicePrivateKey();
    const auto pub_a = instance.alicePublicKey();

    std::cout << "========================================================" << std::endl;
    std::cout << "|alice_private_key| = " << priv_a.length() << std::endl;
    std::cout << "|alice_public_key| = " << pub_a.length() << std::endl;

    // Some diagnostics
    typedef ThLeftNormalForm NF;
    const auto n = protocol.parameters().n();

    if (NF(n, priv_a * -pub_a).isTrivial()) {
      std::cout << "Error: Same Public and private keys!" << std::endl;
      continue;
    }

    std::mt19937_64 g(i);

    int spent_time = time(0);

    const auto is_successful_attack =
        kayawood::attack<decltype(protocol)::field_t, decltype(protocol)::stabilizer_t>(instance, g);

    spent_time = time(0) - spent_time;

    if (is_successful_attack) {
      ++success_count;
    }

    std::cout << "Experiment #" << i << (is_successful_attack ? " is successful." : " failed.") << std::endl;
    std::cout << "Time spent: " << spent_time << std::endl;
    std::cout << "Success: " << success_count << " out of " << (i + 1 - init_seed) << std::endl;
  }

  return 0;
}
