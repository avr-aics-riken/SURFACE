/*
 * LSGL - Large Scale Graphics Library
 *
 * Copyright (c) 2013 Advanced Institute for Computational Science, RIKEN.
 * All rights reserved.
 *
 */

#ifndef __LSGL_ALIGNED_VECTOR_HPP__
#define __LSGL_ALIGNED_VECTOR_HPP__

#include <cassert>
#include <vector>
#include <cstdlib> // posix_memalign

namespace lsgl {
namespace render {

template <typename T, size_t Sz = 16> class aligned_allocator {
public:
#ifdef _WIN32
  static inline void *aligned_malloc(size_t sz, size_t al) {
    return ::_aligned_malloc(sz, al);
  }

  static inline void aligned_free(void *p) { ::_aligned_free(p); }
#elif defined(__unix__) || defined(LINUX) || defined(__APPLE__) ||             \
    defined(__GNUC__)
  static inline void *aligned_malloc(size_t sz, size_t al) {
    void *p = NULL;
    int ret = posix_memalign(&p, al, sz);
    assert(ret == 0);
    return p;
  }
  static inline void aligned_free(void *p) { free(p); }
#endif
public:
  // Typedefs
  typedef size_t size_type;
  typedef ptrdiff_t difference_type;
  typedef T *pointer;
  typedef const T *const_pointer;
  typedef T &reference;
  typedef const T &const_reference;
  typedef T value_type;

public:
  // Constructors
  aligned_allocator() {}
  aligned_allocator(const aligned_allocator &) {}
#if _MSC_VER >= 1300 // VC6 can't handle template members
  template <typename U> aligned_allocator(const aligned_allocator<U> &) {}
#endif
  aligned_allocator &operator=(const aligned_allocator &) { return *this; }

  // Destructor
  ~aligned_allocator() {}

  // Utility functions
  pointer address(reference r) const { return &r; }
  const_pointer address(const_reference c) const { return &c; }
  size_type max_size() const {
    return std::numeric_limits<size_t>::max() / sizeof(T);
  }

  // In-place construction
  void construct(pointer p, const_reference c) {
    new (reinterpret_cast<void *>(p)) T(c); // placement new operator
  }

  // In-place destruction
  void destroy(pointer p) {
    (p)->~T(); // call destructor directly
  }

  // Rebind to allocators of other types
  template <typename U> struct rebind {
    typedef aligned_allocator<U> other;
  };

  // Allocate raw memory
  void deallocate(pointer p, size_type) {
    if (p)
      aligned_free((void *)p);
  }

  pointer allocate(size_type n, const void * = NULL) {
    void *p = aligned_malloc(sizeof(T) * n, Sz);
    assert(p);
    return pointer(p);
  }
}; // end of aligned_allocator

// Comparison
template <typename T1, typename T2>
bool operator==(const aligned_allocator<T1> &, const aligned_allocator<T2> &) {
  return true;
}

template <typename T1, typename T2>
bool operator!=(const aligned_allocator<T1> &, const aligned_allocator<T2> &) {
  return false;
}

namespace aligned_vector_detail {
template <class T> struct same_size_traits {
  typedef struct {
    char buffuer[sizeof(T)];
  } type;
};
}

template <class T, size_t Sz = 16>
class aligned_vector
    : public std::vector<
          typename aligned_vector_detail::same_size_traits<T>::type,
          aligned_allocator<
              typename aligned_vector_detail::same_size_traits<T>::type, Sz>> {
public:
  typedef typename aligned_vector_detail::same_size_traits<T>::type buffer_type;
  typedef std::vector<buffer_type, aligned_allocator<buffer_type>> base_type;

public:
  // aligned_vector(){}
  // aligned_vector(const aligned_vector<T>& rhs):base_type(rhs){}
  T &operator[](size_t i) {
    return reinterpret_cast<T &>(base_type::operator[](i));
  }
  const T &operator[](size_t i) const {
    return reinterpret_cast<const T &>(base_type::operator[](i));
  }
  void push_back(const T &val) {
    base_type::push_back(reinterpret_cast<const buffer_type &>(val));
  }
  T &front() { return reinterpret_cast<T &>(base_type::front()); }
  const T &front() const {
    return reinterpret_cast<const T &>(base_type::front());
  }
  T &back() { return reinterpret_cast<T &>(base_type::back()); }
  const T &back() const {
    return reinterpret_cast<const T &>(base_type::back());
  }
  // void swap(aligned_vector<T>& rhs){base_type::swap((base_type&)rhs);}
};

} // namespace render
} // namespace lsgl
#endif
