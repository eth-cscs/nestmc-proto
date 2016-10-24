/*
 * Author: bcumming
 *
 * Created on May 3, 2014, 5:14 PM
 */

#pragma once

#include <iostream>
#include <type_traits>

#include "definitions.hpp"
#include "util.hpp"
#include "ArrayView.hpp"

////////////////////////////////////////////////////////////////////////////////
namespace memory{

// forward declarations
template<typename T, typename Coord>
class Array;

template <types::size_type Alignment>
class AlignedPolicy;

template<typename T, typename Policy>
class Allocator;

template <typename T, class Allocator>
class HostCoordinator;

namespace util {
    template <typename T, typename Coord>
    struct type_printer<Array<T,Coord>>{
        static std::string print() {
            std::stringstream str;
#if VERBOSE > 1
            str << util::white("Array") << "<" << type_printer<T>::print()
                << ", " << type_printer<Coord>::print() << ">";
#else
            str << util::white("Array") << "<"
                << type_printer<Coord>::print() << ">";
#endif
            return str.str();
        }
    };

    template <typename T, typename Coord>
    struct pretty_printer<Array<T,Coord>>{
        static std::string print(const Array<T,Coord>& val) {
            std::stringstream str;
            str << type_printer<Array<T,Coord>>::print()
                << "(size="     << val.size()
                << ", pointer=" << util::print_pointer(val.data()) << ")";
            return str.str();
        }
    };
}

namespace impl {
    // metafunctions for checking array types
    template <typename T>
    struct is_array_by_value : std::false_type {};

    template <typename T, typename Coord>
    struct is_array_by_value<Array<T, Coord> > : std::true_type {};

    template <typename T>
    struct is_array :
        std::conditional<
            impl::is_array_by_value  <typename std::decay<T>::type> ::value ||
            impl::is_array_view<typename std::decay<T>::type> ::value,
            std::true_type, std::false_type
        >::type
    {};

    // template <typename T, typename Coord>
    // struct is_array<Array<T, Coord> > : std::true_type {};

    template <typename T>
    using is_array_t = typename is_array<T>::type;
}

using impl::is_array;

// array by value
// this wrapper owns the memory in the array
// and is responsible for allocating and freeing memory
template <typename T, typename Coord>
class Array :
    public ArrayView<T, Coord> {
public:
    using value_type = T;
    using base       = ArrayView<value_type, Coord>;
    using view_type  = ArrayView<value_type, Coord>;
    using const_view_type = ConstArrayView<value_type, Coord>;

    using coordinator_type = typename Coord::template rebind<value_type>;

    using typename base::size_type;
    using typename base::difference_type;

    using typename base::pointer;
    using typename base::const_pointer;

    using typename base::iterator;
    using typename base::const_iterator;

    // TODO what about specialized references for things like GPU memory?
    using reference       = value_type&;
    using const_reference = value_type const&;

    // default constructor
    // create empty storage
    Array() :
        base(nullptr, 0)
    {}

    // constructor by size
    template <
        typename I,
        typename = typename std::enable_if<std::is_integral<I>::value>::type>
    Array(I n) :
        base(coordinator_type().allocate(n))
    {
#ifdef VERBOSE
        std::cerr << util::green("Array(" + std::to_string(n) + ")")
                  << "\n  this  " << util::pretty_printer<Array>::print(*this) << std::endl;
#endif
    }

    // constructor by size with default value
    template <
        typename II, typename TT,
        typename = typename std::enable_if<std::is_integral<II>::value>::type,
        typename = typename std::enable_if<std::is_convertible<TT,value_type>::value>::type >
    Array(II n, TT value) :
        base(coordinator_type().allocate(n))
    {
        #ifdef VERBOSE
        std::cerr << util::green("Array(" + std::to_string(n) + ", " + std::to_string(value) + ")")
                  << "\n  this  " << util::pretty_printer<Array>::print(*this) << "\n";
        #endif
        coordinator_type().set(*this, value_type(value));
    }

    // copy constructor
    Array(const Array& other) :
        base(coordinator_type().allocate(other.size()))
    {
        static_assert(impl::is_array_t<Array>::value, "oooooooooooo.............");
#ifdef VERBOSE
        std::cerr << util::green("Array(Array&)")
                  << " " << util::type_printer<Array>::print()
                  << "\n  this  " << util::pretty_printer<Array>::print(*this)
                  << "\n  other " << util::pretty_printer<Array>::print(other) << "\n";
#endif
        coordinator_.copy(const_view_type(other), view_type(*this));
    }

    // move constructor
    Array(Array&& other) {
#ifdef VERBOSE
        std::cerr << util::green("Array(Array&&)")
                  << " " << util::type_printer<Array>::print()
                  << "\n  other " << util::pretty_printer<Array>::print(other) << "\n";
#endif
        base::swap(other);
    }

    // copy constructor where other is an array, array_view or array_reference
    template <
        typename Other,
        typename = typename std::enable_if<impl::is_array_t<Other>::value>::type
    >
    Array(const Other& other) :
        base(coordinator_type().allocate(other.size()))
    {
#ifdef VERBOSE
        std::cerr << util::green("Array(Other&)")
                  << " " << util::type_printer<Array>::print()
                  << "\n  this  " << util::pretty_printer<Array>::print(*this)
                  << "\n  other " << util::pretty_printer<Other>::print(other) << std::endl;
#endif
        coordinator_.copy(typename Other::const_view_type(other), view_type(*this));
    }

    Array& operator=(const Array& other) {
#ifdef VERBOSE
        std::cerr << util::green("Array operator=(Array&)")
                  << "\n  this  "  << util::pretty_printer<Array>::print(*this)
                  << "\n  other " << util::pretty_printer<Array>::print(other) << "\n";
#endif
        coordinator_.free(*this);
        auto ptr = coordinator_type().allocate(other.size());
        base::reset(ptr.data(), other.size());
        coordinator_.copy(const_view_type(other), view_type(*this));
        return *this;
    }

    Array& operator = (Array&& other) {
#ifdef VERBOSE
        std::cerr << util::green("Array operator=(Array&&)")
                  << "\n  this  "  << util::pretty_printer<Array>::print(*this)
                  << "\n  other " << util::pretty_printer<Array>::print(other) << "\n";
#endif
        base::swap(other);
        return *this;
    }

    // have to free the memory in a "by value" range
    ~Array() {
#ifdef VERBOSE
        std::cerr << util::red("~") + util::green("Array()")
                  << "\n  this " << util::pretty_printer<Array>::print(*this) << "\n";
#endif
        coordinator_.free(*this);
    }

    // use the accessors provided by ArrayView
    // this enforces the requirement that accessing all of or a sub-array of an
    // Array should return a view, not a new array.
    using base::operator();

    const coordinator_type& coordinator() const {
        return coordinator_;
    }

    using base::size;

    using base::alignment;

private:

    coordinator_type coordinator_;
};

} // namespace memory
////////////////////////////////////////////////////////////////////////////////

