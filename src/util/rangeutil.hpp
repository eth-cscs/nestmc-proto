#pragma once

/*
 * Sequence and container utilities compatible
 * with ranges.
 */

#include <algorithm>
#include <iterator>
#include <ostream>
#include <numeric>

#include <util/meta.hpp>
#include <util/range.hpp>
#include <util/transform.hpp>
#include <util/meta.hpp>

namespace nest {
namespace mc {
namespace util {

// Present a single item as a range

template <typename T>
range<T*> singleton_view(T& item) {
    return {&item, &item+1};
}

template <typename T>
range<const T*> singleton_view(const T& item) {
    return {&item, &item+1};
}

// Non-owning views and subviews

template <typename Seq>
range<typename sequence_traits<Seq>::iterator, typename sequence_traits<Seq>::sentinel>
range_view(Seq& seq) {
    return make_range(std::begin(seq), std::end(seq));
}

template <
    typename Seq,
    typename Iter = typename sequence_traits<Seq>::iterator,
    typename Size = typename sequence_traits<Seq>::size_type
>
enable_if_t<is_forward_iterator<Iter>::value, range<Iter>>
subrange_view(Seq& seq, Size bi, Size ei) {
    Iter b = std::next(std::begin(seq), bi);
    Iter e = std::next(b, ei-bi);
    return make_range(b, e);
}

// Append sequence to a container

template <typename Container, typename Seq>
Container& append(Container &c, const Seq& seq) {
    auto canon = canonical_view(seq);
    c.insert(c.end(), std::begin(canon), std::end(canon));
    return c;
}

// Assign sequence to a container

template <typename AssignableContainer, typename Seq>
AssignableContainer& assign(AssignableContainer& c, const Seq& seq) {
    auto canon = canonical_view(seq);
    c.assign(std::begin(canon), std::end(canon));
    return c;
}


// Assign sequence to a container with transform `proj`

template <typename AssignableContainer, typename Seq, typename Proj>
AssignableContainer& assign_by(AssignableContainer& c, const Seq& seq, const Proj& proj) {
    assign(c, transform_view(seq, proj));
    return c;
}

// Sort in-place
// Note that a const range reference may wrap non-const iterators.

template <typename Seq>
enable_if_t<!std::is_const<typename sequence_traits<Seq>::reference>::value>
sort(Seq& seq) {
    auto canon = canonical_view(seq);
    std::sort(std::begin(canon), std::end(canon));
}

template <typename Seq>
enable_if_t<!std::is_const<typename sequence_traits<Seq>::reference>::value>
sort(const Seq& seq) {
    auto canon = canonical_view(seq);
    std::sort(std::begin(canon), std::end(canon));
}

// Sort in-place by projection `proj`

template <typename Seq, typename Proj>
enable_if_t<!std::is_const<typename sequence_traits<Seq>::reference>::value>
sort_by(Seq& seq, const Proj& proj) {
    using value_type = typename sequence_traits<Seq>::value_type;
    auto canon = canonical_view(seq);

    std::sort(std::begin(canon), std::end(canon),
        [&proj](const value_type& a, const value_type& b) {
            return proj(a) < proj(b);
        });
}

template <typename Seq, typename Proj>
enable_if_t<!std::is_const<typename sequence_traits<Seq>::reference>::value>
sort_by(const Seq& seq, const Proj& proj) {
    using value_type = typename sequence_traits<Seq>::value_type;
    auto canon = canonical_view(seq);

    std::sort(std::begin(canon), std::end(canon),
        [&proj](const value_type& a, const value_type& b) {
            return proj(a) < proj(b);
        });
}

// Stable sort in-place by projection `proj`

template <typename Seq, typename Proj>
enable_if_t<!std::is_const<typename sequence_traits<Seq>::reference>::value>
stable_sort_by(Seq& seq, const Proj& proj) {
    using value_type = typename sequence_traits<Seq>::value_type;
    auto canon = canonical_view(seq);

    std::stable_sort(std::begin(canon), std::end(canon),
        [&proj](const value_type& a, const value_type& b) {
            return proj(a) < proj(b);
        });
}

template <typename Seq, typename Proj>
enable_if_t<!std::is_const<typename sequence_traits<Seq>::reference>::value>
stable_sort_by(const Seq& seq, const Proj& proj) {
    using value_type = typename sequence_traits<Seq>::value_type;
    auto canon = canonical_view(seq);

    std::stable_sort(std::begin(canon), std::end(canon),
        [&proj](const value_type& a, const value_type& b) {
            return proj(a) < proj(b);
        });
}

// Range-interface for `all_of`, `any_of`

template <typename Seq, typename Predicate>
bool all_of(const Seq& seq, const Predicate& pred) {
    auto canon = canonical_view(seq);
    return std::all_of(std::begin(canon), std::end(canon), pred);
}

template <typename Seq, typename Predicate>
bool any_of(const Seq& seq, const Predicate& pred) {
    auto canon = canonical_view(seq);
    return std::any_of(std::begin(canon), std::end(canon), pred);
}


// Accumulate by projection `proj`

template <
    typename Seq,
    typename Proj,
    typename Value = typename transform_iterator<typename sequence_traits<Seq>::const_iterator, Proj>::value_type
>
Value sum_by(const Seq& seq, const Proj& proj, Value base = Value{}) {
    auto canon = canonical_view(transform_view(seq, proj));
    return std::accumulate(std::begin(canon), std::end(canon), base);
}

// Maximum element by projection `proj`
// - returns an iterator `i` into supplied sequence which has the maximum
//   value of `proj(*i)`.

template <typename Seq, typename Proj>
typename sequence_traits<Seq>::iterator
max_element_by(Seq& seq, const Proj& proj) {
    using value_type = typename sequence_traits<Seq>::value_type;
    auto canon = canonical_view(seq);

    return std::max_element(std::begin(canon), std::end(canon),
        [&proj](const value_type& a, const value_type& b) {
            return proj(a) < proj(b);
        });
}

template <typename Seq, typename Proj>
typename sequence_traits<Seq>::iterator
max_element_by(const Seq& seq, const Proj& proj) {
    using value_type = typename sequence_traits<Seq>::value_type;
    auto canon = canonical_view(seq);

    return std::max_element(std::begin(canon), std::end(canon),
        [&proj](const value_type& a, const value_type& b) {
            return proj(a) < proj(b);
        });
}

// Maximum value
//
// Value semantics instead of iterator semantics means it will operate
// with input iterators.  Will return default-constructed value if sequence
// is empty.
//
// (Consider making generic associative reduction with TBB implementation
// for random-access iterators?)

template <
    typename Seq,
    typename Value = typename sequence_traits<Seq>::value_type,
    typename Compare = std::less<Value>
>
Value max_value(const Seq& seq, Compare cmp = Compare{}) {
    if (util::empty(seq)) {
        return Value{};
    }

    auto i = std::begin(seq);
    auto e = std::end(seq);
    auto m = *i;
    while (++i!=e) {
        Value x = *i;
        if (cmp(m, x)) {
            m = std::move(x);
        }
    }
    return m;
}

template <typename T, typename Seq>
std::vector<T> make_std_vector(const Seq& seq) {
    auto i = std::begin(seq);
    auto e = std::end(seq);
    return std::vector<T>(i, e);
}

template <typename C, typename Seq>
C make_copy(Seq const& seq) {
    return C{std::begin(seq), std::end(seq)};
}

} // namespace util
} // namespace mc
} // namespace nest

