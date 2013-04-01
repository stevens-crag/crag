// Copyright 2006, Google Inc.
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
//     * Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above
// copyright notice, this list of conditions and the following disclaimer
// in the documentation and/or other materials provided with the
// distribution.
//     * Neither the name of Google Inc. nor the names of its
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <iostream>
#include <cstring>
#include <gmp.h>

#include <boost/pool/singleton_pool.hpp>
//
#include "gtest/gtest.h"

//namespace {

struct GmpPoolTag {};
const size_t POOL_ALLOCATE_BOUNDARY = 4*sizeof(mp_limb_t);
size_t allocate_pool = 0;
size_t allocate_malloc = 0;
size_t reallocate_pool = 0;
size_t reallocate_pool_malloc = 0;
size_t reallocate_malloc = 0;
size_t free_pool = 0;
size_t free_malloc = 0;

void* (*gmp_default_malloc) (size_t);
void* (*gmp_default_realloc) (void *, size_t, size_t);
void (*gmp_default_free) (void *, size_t);

typedef boost::singleton_pool<GmpPoolTag, POOL_ALLOCATE_BOUNDARY, boost::default_user_allocator_malloc_free, boost::details::pool::null_mutex> GmpPool;

inline void* gmp_pool_allocate(size_t alloc_size) {
  if (alloc_size == 0) {
    return gmp_default_malloc(alloc_size);
  } else if (alloc_size <= POOL_ALLOCATE_BOUNDARY) {
    //++allocate_pool;
    return GmpPool::malloc();
  } else {
    //++allocate_malloc;
    return gmp_default_malloc(alloc_size);
  }
}

inline void* gmp_pool_reallocate(void *ptr, size_t old_size, size_t new_size) {
  if (new_size == 0) {
    return nullptr;
  } else if (new_size <= POOL_ALLOCATE_BOUNDARY) {
    if (old_size) {
      return ptr;
    } else {
      //++reallocate_pool;
      return GmpPool::malloc();
    }
  } else if (old_size == 0) {
    //++reallocate_malloc;
    return gmp_default_realloc(ptr, old_size, new_size);
  } else if (old_size <= POOL_ALLOCATE_BOUNDARY) {
    //++reallocate_pool_malloc;
    void * new_block = malloc(new_size);
    memmove(new_block, ptr, old_size);
    GmpPool::free(ptr);
    return new_block;
  } else {
    //++reallocate_malloc;
    return gmp_default_realloc(ptr, old_size, new_size);
  }
}

inline void gmp_pool_free(void *ptr, size_t size) {
  if (size <= POOL_ALLOCATE_BOUNDARY) {
    //++free_pool;
    GmpPool::free(ptr);
  } else {
    //++free_malloc;
    gmp_default_free(ptr, size);
  }
}

//} //anonymous namespace
GTEST_API_ int main(int argc, char **argv) {
  mp_get_memory_functions (&gmp_default_malloc, &gmp_default_realloc, &gmp_default_free);
  mp_set_memory_functions(gmp_pool_allocate, gmp_pool_reallocate, gmp_pool_free);
  //mp_set_memory_functions(gmp_default_malloc, gmp_default_realloc, gmp_default_free);
  std::cout << "Running main() from gtest_main.cc\n";

  testing::InitGoogleTest(&argc, argv);
  bool ret = RUN_ALL_TESTS();
  std::cout << allocate_pool << std::endl;
  std::cout << allocate_malloc << std::endl;
  std::cout << reallocate_pool << std::endl;
  std::cout << reallocate_pool_malloc << std::endl;
  std::cout << reallocate_malloc << std::endl;
  std::cout << free_pool << std::endl;
  std::cout << free_malloc << std::endl;
  return ret;
}
