/*
 * gmp_boost_pool_allocator.h
 *
 *  Created on: Apr 2, 2013
 *      Author: dpantele
 */

#ifndef GMP_BOOST_POOL_ALLOCATOR_H_
#define GMP_BOOST_POOL_ALLOCATOR_H_

#include <cstring>
#include <gmp.h>
#include "boost/pool/singleton_pool.hpp"

namespace {

struct GmpPoolTag {};
const size_t POOL_ALLOCATE_BOUNDARY = 4 * sizeof(mp_limb_t);

extern void* (*gmp_default_malloc) (size_t);
extern void* (*gmp_default_realloc) (void *, size_t, size_t);
extern void (*gmp_default_free) (void *, size_t);

typedef boost::singleton_pool<GmpPoolTag, POOL_ALLOCATE_BOUNDARY, boost::default_user_allocator_malloc_free, boost::details::pool::null_mutex> GmpPool;

inline void* gmp_pool_allocate(size_t alloc_size) {
  if (alloc_size <= POOL_ALLOCATE_BOUNDARY) {
    return GmpPool::malloc();
  } else {
    return gmp_default_malloc(alloc_size);
  }
}

inline void* gmp_pool_reallocate(void *ptr, size_t old_size, size_t new_size) {
  if (new_size <= POOL_ALLOCATE_BOUNDARY) {
    if (old_size) {
      return ptr;
    } else {
      return GmpPool::malloc();
    }
  } else if (old_size == 0) {
    return gmp_default_realloc(ptr, old_size, new_size);
  } else if (old_size <= POOL_ALLOCATE_BOUNDARY) {
    void * new_block = malloc(new_size);
    memmove(new_block, ptr, old_size);
    GmpPool::free(ptr);
    return new_block;
  } else {
    return gmp_default_realloc(ptr, old_size, new_size);
  }
}

inline void gmp_pool_free(void *ptr, size_t size) {
  if (size <= POOL_ALLOCATE_BOUNDARY) {
    GmpPool::free(ptr);
  } else {
    gmp_default_free(ptr, size);
  }
}

} //end of local namespace

void gmp_pool_setup();

#endif /* GMP_BOOST_POOL_ALLOCATOR_H_ */
