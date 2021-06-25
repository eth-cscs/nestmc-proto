#pragma once

/* Present a pair of iterators as a non-owning collection.
 *
 * Two public member fields, `left` and `right`, describe
 * the half-open interval [`left`, `right`).
 *
 * Constness of the range object only affects mutability
 * of the iterators, and does not relate to the constness
 * of the data to which the iterators refer.
 *
 * The `right` field may differ in type from the `left` field,
 * in which case it is regarded as a sentinel type; the end of
 * the interval is then marked by the first successor `i` of
 * `left` that satisfies `i==right`.
 *
 * For an iterator `i` and sentinel `s`, it is expected that
 * the tests `i==s` and `i!=s` are well defined, with the
 * corresponding semantics.
 */

#include <cstddef>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>

#include <arbor/assert.hpp>

#include "util/counter.hpp"
#include "util/iterutil.hpp"
#include "util/meta.hpp"
#include "util/sentinel.hpp"

namespace arb {
namespace util {

template <typename U, typename S = U>
struct range {
    using iterator = U;
    using sentinel = S;
    using const_iterator = iterator;
    using difference_type = typename std::iterator_traits<iterator>::difference_type;
    using size_type = std::make_unsigned_t<difference_type>;
    using value_type = typename std::iterator_traits<iterator>::value_type;
    using reference = typename std::iterator_traits<iterator>::reference;
    using const_reference = const value_type&;

    iterator left;
    sentinel right;

    range() = default;
    range(const range&) = default;
    range(range&&) = default;

    template <typename U1, typename U2>
    range(U1&& l, U2&& r):
        left(std::forward<U1>(l)), right(std::forward<U2>(r))
    {}

    template <typename U1, typename U2>
    range(const std::pair<U1, U2>& p):
        left(p.first), right(p.second)
    {}

    template <
        typename U1,
        typename U2,
        typename = std::enable_if_t<
            std::is_constructible<iterator, U1>::value &&
            std::is_constructible<sentinel, U2>::value>
    >
    range(const range<U1, U2>& other):
        left(other.left), right(other.right)
    {}

    range& operator=(const range&) = default;
    range& operator=(range&&) = default;

    template <typename U1, typename U2>
    range& operator=(const range<U1, U2>& other) {
        left = other.left;
        right = other.right;
        return *this;
    }

    bool empty() const { return left == right; }

    iterator begin() const { return left; }
    const_iterator cbegin() const { return left; }

    sentinel end() const { return right; }
    sentinel cend() const { return right; }

    template <typename V = iterator>
    std::enable_if_t<is_forward_iterator<V>::value, size_type>
    size() const {
        return util::distance(begin(), end());
    }

    constexpr size_type max_size() const {
        return std::numeric_limits<size_type>::max();
    }

    void swap(range& other) {
        std::swap(left, other.left);
        std::swap(right, other.right);
    }

    decltype(auto) front() const { return *left; }

    decltype(auto) back() const { return *upto(left, right); }

    template <typename V = iterator>
    std::enable_if_t<is_random_access_iterator<V>::value, decltype(*left)>
    operator[](difference_type n) const {
        return *std::next(begin(), n);
    }

    template <typename V = iterator>
    std::enable_if_t<is_random_access_iterator<V>::value, decltype(*left)>
    at(difference_type n) const {
        if (size_type(n) >= size()) {
            throw std::out_of_range("out of range in range");
        }
        return (*this)[n];
    }

    // Expose `data` method if a pointer range.
    template <typename V = iterator, typename W = sentinel>
    std::enable_if_t<std::is_same<V, W>::value && std::is_pointer<V>::value, iterator>
    data() const {
        return left;
    }
};

template <typename U, typename V>
range<U, V> make_range(const U& left, const V& right) {
    return range<U, V>(left, right);
}

template <typename U, typename V>
range<U, V> make_range(const std::pair<U, V>& iterators) {
    return range<U, V>(iterators.first, iterators.second);
}

// From pointer and length
template <typename T>
range<T, T> range_n(T t, size_t n) {
    return {t, t + n};
}

// Present a possibly sentinel-terminated range as an STL-compatible sequence
// using the sentinel_iterator adaptor.

template <typename Seq>
auto canonical_view(Seq&& s) {
    using std::begin;
    using std::end;

    return make_range(
        make_sentinel_iterator(begin(s), end(s)),
        make_sentinel_end(begin(s), end(s)));
}

// Strictly evaluate end point in sentinel-terminated range and present as a range over
// iterators. Note: O(N) behaviour with forward iterator ranges or sentinel-terminated ranges.

template <typename Seq>
auto strict_view(Seq&& s) {
    using std::begin;
    using std::end;

    auto b = begin(s);
    auto e = end(s);
    return make_range(b, b==e? b: std::next(util::upto(b, e)));
}

} // namespace util
} // namespace arb
