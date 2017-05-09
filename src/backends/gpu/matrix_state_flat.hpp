#pragma once

#include <memory/memory.hpp>
#include <util/span.hpp>
#include <util/partition.hpp>
#include <util/rangeutil.hpp>

#include "kernels/solve_matrix.hpp"
#include "kernels/assemble_matrix.hpp"

namespace nest {
namespace mc {
namespace gpu {

/// matrix state
template <typename T, typename I>
struct matrix_state_flat {
    using value_type = T;
    using size_type = I;

    using array  = memory::device_vector<value_type>;
    using iarray = memory::device_vector<size_type>;

    using view = typename array::view_type;
    using const_view = typename array::const_view_type;

    iarray parent_index;
    iarray cell_cv_divs;
    iarray cv_to_cell;

    array d;     // [μS]
    array u;     // [μS]
    array rhs;   // [nA]

    array cv_capacitance;      // [pF]
    array face_conductance;    // [μS]

    // the invariant part of the matrix diagonal
    array invariant_d;         // [μS]

    // interface for exposing the solution to the outside world
    view solution;

    matrix_state_flat() = default;

    matrix_state_flat(const std::vector<size_type>& p,
                 const std::vector<size_type>& cell_cv_divs,
                 const std::vector<value_type>& cv_cap,
                 const std::vector<value_type>& face_cond):
        parent_index(memory::make_const_view(p)),
        cell_cv_divs(memory::make_const_view(cell_cv_divs)),
        cv_to_cell(p.size()),
        d(p.size()),
        u(p.size()),
        rhs(p.size()),
        cv_capacitance(memory::make_const_view(cv_cap))
    {
        EXPECTS(cv_cap.size() == size());
        EXPECTS(face_cond.size() == size());
        EXPECTS(cell_cv_divs.back() == size());
        EXPECTS(cell_cv_divs.size() > 2u);

        using memory::make_const_view;

        auto n = d.size();
        std::vector<size_type> cv_to_cell_tmp(n, 0);
        std::vector<value_type> invariant_d_tmp(n, 0);
        std::vector<value_type> u_tmp(n, 0);

        for (auto i: util::make_span(1u, n)) {
            auto gij = face_cond[i];

            u_tmp[i] = -gij;
            invariant_d_tmp[i] += gij;
            invariant_d_tmp[p[i]] += gij;
        }

        size_type ci = 0;
        for (auto cv_span: util::partition_view(cell_cv_divs)) {
            util::fill(util::subrange_view(cv_to_cell_tmp, cv_span), ci);
            ++ci;
        }

        cv_to_cell = make_const_view(cv_to_cell_tmp);
        invariant_d = make_const_view(invariant_d_tmp);
        u = make_const_view(u_tmp);

        solution = rhs;
    }

    // Assemble the matrix
    // Afterwards the diagonal and RHS will have been set given dt, voltage and current,
    // where dt is determined by the start and end integration times t and t_to.
    //   t       [ms]
    //   t_to    [ms]
    //   voltage [mV]
    //   current [nA]
    void assemble(const_view t, const_view t_to, const_view voltage, const_view current) {
        // determine the grid dimensions for the kernel
        auto const n = voltage.size();
        auto const block_dim = 128;
        auto const grid_dim = impl::block_count(n, block_dim);

        assemble_matrix_flat<value_type, size_type><<<grid_dim, block_dim>>> (
            d.data(), rhs.data(), invariant_d.data(), voltage.data(),
            current.data(), cv_capacitance.data(), cv_to_cell.data(), t.data(), t_to.data(), size());
    }

    void solve() {
        // determine the grid dimensions for the kernel
        auto const block_dim = 128;
        auto const grid_dim = impl::block_count(num_matrices(), block_dim);

        // perform solve on gpu
        solve_matrix_flat<value_type, size_type><<<grid_dim, block_dim>>> (
            rhs.data(), d.data(), u.data(), parent_index.data(),
            cell_cv_divs.data(), num_matrices());
    }

    std::size_t size() const {
        return parent_index.size();
    }

private:
    unsigned num_matrices() const {
        return cell_cv_divs.size()-1;
    }
};

} // namespace gpu
} // namespace mc
} // namespace nest
