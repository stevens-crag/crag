#include <gtest/gtest.h>

#include "kayawood_attack.h"

namespace crag {
namespace kayawood {
namespace {

TEST(KayawoodAttack, Params128bit) {
  const auto protocol = kayawood::getProtocolFor128BitsSecurity(0);
  const auto instance = protocol.generateInstance(0);
}

TEST(KayawoodAttack, Params128bitMultipleCloaking) {
  const auto protocol = kayawood::getProtocolFor128BitsSecurityMultipleCloaking(0);
  const auto instance = protocol.generateInstance(0);
}

TEST(KayawoodAttack, Params256bit) {
  const auto protocol = kayawood::getProtocolFor256BitsSecurity(0);
  const auto instance = protocol.generateInstance(0);
}

TEST(KayawoodAttack, Params256bitMultipleCloaking) {
  const auto protocol = kayawood::getProtocolFor256BitsSecurityMultipleCloaking(0);
  const auto instance = protocol.generateInstance(0);
}
} // namespace
} // namespace kayawood
} // namespace crag
