#pragma once

#ifndef CRAG_PARALLEL_H
#define CRAG_PARALLEL_H

#include <boost/container/vector.hpp>

#include <atomic>
#include <thread>
#include <vector>

namespace crag {
namespace parallel {

size_t getHardwareConcurrency();

//! Parallel for-each invoking fn for each i = 0,...,n-1.
//! Function is of type
//!     void fn(size_t)
template <typename Function>
void forEach(size_t n, Function fn) {
  std::atomic<size_t> index(0);

  const auto thread_fn = [&]() {
    while (true) {
      const auto i = index.fetch_add(1);

      if (i >= n) {
        return;
      }

      fn(i);
    }
  };

  const auto threads_count = std::min(n, getHardwareConcurrency());

  std::vector<std::thread> threads;
  threads.reserve(threads_count);

  for (size_t i = 0; i < threads_count; ++i) {
    threads.push_back(std::thread(thread_fn));
  }

  for (size_t i = 0; i < threads_count; ++i) {
    threads[i].join();
  }
}

//! Parallel for-each.
//! Function is of type
//!     void fn(size_t, const T&) or
//!     void fn(size_t, T)
template <typename T, typename Function>
void forEach(const std::vector<T>& items, Function fn) {
  const auto size = items.size();

  forEach(size, [&](size_t i) { fn(i, items[i]); });
}

//! Parallel map, TOut must be default constructible.
//! Function is of type
//!     TOut fn(const TIn&) or
//!     TOut fn(TIn)
template <
    typename TIn,
    typename TOut,
    typename Function,
    typename = typename std::enable_if<!std::is_same<TOut, bool>::value>::type>
std::vector<TOut> map(const std::vector<TIn>& items, Function fn) {
  std::vector<TOut> result(items.size());

  forEach(items, [&](size_t i, const TIn& item) { result[i] = fn(item); });

  return result;
}

//! Workaround for std::vector<bool>
template <typename T, typename Function>
boost::container::vector<bool> bmap(const std::vector<T>& items, Function fn) {
  boost::container::vector<bool> result(items.size());

  forEach(items.size(), [&](size_t i) { result[i] = fn(items[i]); });

  return result;
}

//! Workaround for std::vector<bool>
template <typename Function>
boost::container::vector<bool> bmap(size_t n, Function fn) {
  boost::container::vector<bool> result(n);

  forEach(n, [&](size_t i) { result[i] = fn(i); });

  return result;
}

//! Parallel map, T must be default constructible.
//! Function is of type
//!     T fn(const T&) or
//!     T fn(T)
template <typename T, typename Function, typename = typename std::enable_if<!std::is_same<T, bool>::value>::type>
std::vector<T> map(const std::vector<T>& items, Function fn) {
  return map<T, T>(items, std::move(fn));
}

//! Parallel map, T must be default constructible.
//! Function is of type
//!     T fn(size_t)
template <typename T, typename Function, typename = typename std::enable_if<!std::is_same<T, bool>::value>::type>
std::vector<T> map(size_t n, Function fn) {
  std::vector<T> result(n);

  forEach(n, [&](size_t i) { result[i] = fn(i); });

  return result;
}
} // namespace parallel
} // namespace crag

#endif // CRAG_PARALLEL_H
