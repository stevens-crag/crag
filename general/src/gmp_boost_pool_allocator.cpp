#include "gmp_boost_pool_allocator.h"

namespace {

void* (*gmp_default_malloc) (size_t) = nullptr;
void* (*gmp_default_realloc) (void *, size_t, size_t) = nullptr;
void (*gmp_default_free) (void *, size_t) = nullptr;

}

void gmp_pool_setup() {
  mp_get_memory_functions (&gmp_default_malloc, &gmp_default_realloc, &gmp_default_free);
  mp_set_memory_functions(gmp_pool_allocate, gmp_pool_reallocate, gmp_pool_free);
}
