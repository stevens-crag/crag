#include "gmp_boost_pool_allocator.h"

#include <cstring>
#include "gmp.h"

#include "boost/pool/singleton_pool.hpp"

namespace gmp_boost_pool_detail {

struct GmpPoolTag {};
const size_t POOL_ALLOCATE_BOUNDARY = 4 * sizeof(mp_limb_t);

typedef boost::singleton_pool<
  GmpPoolTag,
  POOL_ALLOCATE_BOUNDARY,
  boost::default_user_allocator_malloc_free,
  boost::details::pool::null_mutex
> GmpPool;

void* (*gmp_default_malloc) (size_t) = nullptr;
void* (*gmp_default_realloc) (void *, size_t, size_t) = nullptr;
void (*gmp_default_free) (void *, size_t) = nullptr;

inline void* allocate(size_t alloc_size) {
  if (alloc_size <= POOL_ALLOCATE_BOUNDARY) {
    return GmpPool::malloc();
  } else {
    return gmp_default_malloc(alloc_size);
  }
}

inline void* reallocate(void *ptr, size_t old_size, size_t new_size) {
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

inline void free(void *ptr, size_t size) {
  if (size <= POOL_ALLOCATE_BOUNDARY) {
    GmpPool::free(ptr);
  } else {
    gmp_default_free(ptr, size);
  }
}

} //namespace gmp_boost_pool_detail

void gmp_pool_setup() {
  mp_get_memory_functions(
    &gmp_boost_pool_detail::gmp_default_malloc,
    &gmp_boost_pool_detail::gmp_default_realloc,
    &gmp_boost_pool_detail::gmp_default_free
  );

  mp_set_memory_functions(
    gmp_boost_pool_detail::allocate,
    gmp_boost_pool_detail::reallocate,
    gmp_boost_pool_detail::free
  );
}
