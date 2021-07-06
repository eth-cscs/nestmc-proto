#pragma once

// Base class for scoped-instrumentation of glibc malloc.
//
// For the lifetime of a `with_instrumented_malloc` object,
// global memory allocation hooks will be set so that
// the virtual `on_malloc`, `on_realloc`, `on_memalign`
// and `on_free` calls will be invoked before the corresponding
// `malloc`, `realloc` etc. is executed.
//
// Scopes of `with_instrumented_malloc` may be nested, but:
//   * Don't interleave lifetimes of these objects and expect things
//     to work!
//   * Don't try and create new `with_instrumented_malloc` instances
//     from within an `on_malloc` callback (or others).
//   * Definitely don't try and use this in a multithreaded context.
//
// Calling code should check CAN_INSTRUMENT_MALLOC preprocessor
// symbol to see if this functionality is available.

#include <cstddef>

#if (__GLIBC__==2)
#include <malloc.h>
#define CAN_INSTRUMENT_MALLOC
#endif

// Disable if using address sanitizer though:

// This is how clang tells us.
#if defined(__has_feature)
#if __has_feature(address_sanitizer)
#undef CAN_INSTRUMENT_MALLOC
#endif
#  if __has_feature(thread_sanitizer)
#undef CAN_INSTRUMENT_MALLOC
#  endif
#endif
// This is how gcc tells us.
#if defined(__SANITIZE_ADDRESS__)
#undef CAN_INSTRUMENT_MALLOC
#endif

namespace testing {

#ifdef CAN_INSTRUMENT_MALLOC

// For run-time, temporary intervention in the malloc-family calls,
// there is still no better alternative than to use the
// deprecated __malloc_hook pointers and friends. 

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
#if defined(__INTEL_COMPILER)
#pragma warning push
#pragma warning disable 1478
#endif

// Totally not thread safe!
struct with_instrumented_malloc {
    with_instrumented_malloc() {
        push();
    }

    ~with_instrumented_malloc() {
        pop();
    }

    virtual void on_malloc(std::size_t, const void*) {}
    virtual void on_realloc(void*, std::size_t, const void*) {}
    virtual void on_free(void*, const void*) {}
    virtual void on_memalign(std::size_t, std::size_t, const void*) {}

private:
    static with_instrumented_malloc*& instance() {
        static with_instrumented_malloc* ptr = nullptr;
        return ptr;
    }

    with_instrumented_malloc* prev_;
    decltype(__malloc_hook) saved_malloc_hook_;
    decltype(__realloc_hook) saved_realloc_hook_;
    decltype(__free_hook) saved_free_hook_;
    decltype(__memalign_hook) saved_memalign_hook_;

    void push() {
        saved_malloc_hook_ = __malloc_hook;
        saved_realloc_hook_ = __realloc_hook;
        saved_free_hook_ = __free_hook;
        saved_memalign_hook_ = __memalign_hook;

        prev_ = instance();
        instance() = this;

        __malloc_hook = malloc_hook;
        __realloc_hook = realloc_hook;
        __free_hook = free_hook;
        __memalign_hook = memalign_hook;
    }

    void pop() {
        instance() = prev_;
        __malloc_hook = saved_malloc_hook_;
        __realloc_hook = saved_realloc_hook_;
        __free_hook = saved_free_hook_;
        __memalign_hook = saved_memalign_hook_;
    }

    struct windback_guard {
        with_instrumented_malloc* p;

        windback_guard(): p(instance()) { p->pop(); }
        ~windback_guard() { p->push(); }
    };

    static void* malloc_hook(std::size_t size, const void* caller) {
        windback_guard g;
        g.p->on_malloc(size, caller);
        return malloc(size);
    }

    static void* realloc_hook(void* ptr, std::size_t size, const void* caller) {
        windback_guard g;
        g.p->on_realloc(ptr, size, caller);
        return realloc(ptr, size);
    }

    static void free_hook(void* ptr, const void* caller) {
        windback_guard g;
        g.p->on_free(ptr, caller);
        free(ptr);
    }

    static void* memalign_hook(std::size_t alignment, std::size_t size, const void* caller) {
        windback_guard g;
        g.p->on_memalign(alignment, size, caller);
        return memalign(alignment, size);
    }
};

#pragma GCC diagnostic pop
#if defined(__INTEL_COMPILER)
#pragma warning pop
#endif

#else

struct with_instrumented_malloc {
    with_instrumented_malloc() {
        throw std::runtime_error("malloc instrumentation not supported\n");
    }

    virtual void on_malloc(std::size_t, const void*) {}
    virtual void on_realloc(void*, std::size_t, const void*) {}
    virtual void on_free(void*, const void*) {}
    virtual void on_memalign(std::size_t, std::size_t, const void*) {}
};

#endif // ifdef CAN_INSTRUMENT_MALLOC

} // namespace testing
