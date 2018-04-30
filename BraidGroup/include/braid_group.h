#pragma once

#ifndef CRAG_BRAID_GROUP_H
#define CRAG_BRAID_GROUP_H

#include "Word.h"

namespace crag {
namespace braidgroup {

//! Defines a representation of a Braid Group.
//! A BraidGroup is uniquely defined by its index.
class BraidGroup {
public:
  BraidGroup(size_t n)
      : n_(n) {}

  //! get the rank of the braid group (number of strands)
  size_t getRank() const {
    return n_;
  }

  //! get the braid word called a half-twist (its square generates the center of the braid group)
  Word twist(const Word& w) const;

private:
  size_t n_;
};

//! Returns relations in B_n
std::vector<Word> getRelations(size_t n);

//! Given B_n, returns a pure braid subgroup generator for the pair (i, j), where 1 <= i < j <= n.
//! The result is g_{i, j} = b_{j-1} b_{j-2} ... b_{i+1} b_i^2 b_{i+1}^{-1} ... b_{j-2}^{-1} b_{j-1}^{-1}.
Word getPureBraidSubgroupGenerator(size_t n, size_t i, size_t j);
} // namespace braidgroup
} // namespace crag

#endif // CRAG_BRAID_GROUP_H
