#pragma once

// Embedding of cell morphology as 1-d tree with piecewise linear radius.

#include <vector>

#include <arbor/morph/morphology.hpp>
#include <arbor/morph/primitives.hpp>

namespace arb {

struct embed_pwlin_data;

namespace util {
template <typename X> struct pw_elements;
}

// Piecewise-constant functions are represented as scalar (double)
// values defined over contiguous intervals.
using pw_constant_fn = util::pw_elements<double>;

struct embed_pwlin {
    explicit embed_pwlin(const arb::morphology& m);

    // Locations that mark the boundaries between segments.
    const std::vector<mlocation>& segment_locations() const {
        return segment_locations_;
    }

    std::pair<const mlocation*, const mlocation*> branch_segment_locations(msize_t i) const {
        const mlocation* p = segment_locations_.data();
        return std::make_pair(p+branch_segment_part_[i], p+branch_segment_part_[i+1]);
    }

    // Interpolated radius at location.
    double radius(mlocation) const;
    mcable_list radius_cmp(msize_t bid, double rad_lim, comp_op op) const;

    double directed_projection(mlocation) const;
    mcable_list projection_cmp(msize_t bid, double proj_lim, comp_op op) const;

    // Computed length of mcable.
    double integrate_length(mcable c) const;
    double integrate_length(mlocation proxmal, mlocation distal) const;

    double integrate_length(mcable c, const pw_constant_fn&) const;
    double integrate_length(msize_t bid, const pw_constant_fn&) const;

    // Membrane surface area of given mcable.
    double integrate_area(mcable c) const;
    double integrate_area(mlocation proxmal, mlocation distal) const;

    double integrate_area(mcable c, const pw_constant_fn&) const;
    double integrate_area(msize_t bid, const pw_constant_fn&) const;

    // Integrated inverse cross-sectional area of given mcable.
    double integrate_ixa(mcable c) const;

    double integrate_ixa(mcable c, const pw_constant_fn&) const;
    double integrate_ixa(msize_t bid, const pw_constant_fn&) const;

    // Length of whole branch.
    double branch_length(msize_t bid) const {
        return integrate_length(mcable{bid, 0, 1});
    }

private:
    std::vector<mlocation> segment_locations_;
    std::vector<msize_t> branch_segment_part_;
    std::shared_ptr<embed_pwlin_data> data_;
};

} // namespace arb


