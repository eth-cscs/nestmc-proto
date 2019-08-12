#include <algorithm>
#include <stdexcept>
#include <vector>

#include "symge.hpp"

namespace symge {

struct pivot {
    unsigned row;
    unsigned col;
};


// Returns q[c]*p - p[c]*q; new symbols required due to fill-in are provided by the
// `define_sym` functor, which takes a `symbol_term_diff` and returns a `symbol`.

template <typename DefineSym>
sym_row row_reduce(unsigned c, const sym_row& p, const sym_row& q, DefineSym define_sym) {
    std::cout << p.index(c) << " and " << q.index(c) << std::endl;
    if (p.index(c)==p.npos || q.index(c)==q.npos) throw std::runtime_error("improper row reduction");

    sym_row u;
    symbol x = q[c];
    symbol y = p[c];

    auto piter = p.begin();
    auto qiter = q.begin();
    unsigned pj = piter->col;
    unsigned qj = qiter->col;

    while (piter!=p.end() || qiter!=q.end()) {
        unsigned j = std::min(pj, qj);
        symbol_term t1, t2;

        if (j==pj) {
            t1 = x*piter->value;
            ++piter;
            pj = piter==p.end()? p.npos: piter->col;
        }
        if (j==qj) {
            t2 = y*qiter->value;
            ++qiter;
            qj = qiter==q.end()? q.npos: qiter->col;
        }
        if (j!=c) {
            u.push_back({j, define_sym(t1-t2)});
        }
    }
    return u;
}

// Estimate cost of a choice of pivot for G–J reduction below. Uses a simple greedy
// estimate based on immediate fill cost.
double estimate_cost(const sym_matrix& A, pivot p) {
    unsigned nfill = 0;

    auto count_fill = [&nfill](symbol_term_diff t) {
        bool l = t.left;
        bool r = t.right;
        nfill += r&!l;
        return symbol{};
    };

    for (unsigned i = 0; i<A.nrow(); ++i) {
        if (i==p.row || A[i].index(p.col)==msparse::row_npos) continue;
        row_reduce(p.col, A[i], A[p.row], count_fill);
    }

    return nfill;
}

// Perform Gauss-Jordan elimination on given symbolic matrix. New symbols
// required due to fill-in are added to the supplied symbol table.
//
// The matrix A is regarded as being diagonally dominant, and so pivots
// are selected from the diagonal. The choice of pivot at each stage of
// the reduction is goverened by a cost estimation (see above).
//
// The reduction is division-free: the result will have non-zero terms
// that are symbols that are either primitive, or defined (in the symbol
// table) as products or differences of products of other symbols.
void gj_reduce(sym_matrix& A, symbol_table& table) {
    if (A.nrow()>A.ncol()) throw std::runtime_error("improper matrix for reduction");

    auto define_sym = [&table](symbol_term_diff t) { return table.define(t); };

    std::vector<pivot> pivots;
    for (unsigned r = 0; r<A.nrow(); ++r) {
        pivot p;
        p.row = r;
        const sym_row& row = A[r];
        for (unsigned c = 0; c < A.ncol(); ++c) {
            if (row[c] != symbol{}) {
                bool used = false;
                for (auto p: pivots) {
                    used = (p.col == c);
                    if (used) break;
                }
                if (!used) {
                    p.col = c;
                    break;
                }
            }
        }
        std::cout << "(" << p.row << ", " << p.col << ")" << std::endl;
        pivots.push_back(std::move(p));
    }

    for (auto p: pivots) {
        std::cout << "(" << p.row << ", " << p.col << ")" << std::endl;
    }

    std::vector<double> cost(pivots.size());

    while (!pivots.empty()) {
        for (unsigned i = 0; i<pivots.size(); ++i) {
            cost[pivots[i].row] = estimate_cost(A, pivots[i]);
        }

        std::sort(pivots.begin(), pivots.end(),
            [&](pivot r1, pivot r2) { return cost[r1.row]>cost[r2.row]; });

        std::cout << std::endl;
        for (auto p: pivots) {
            std::cout << "(" << p.row << ", " << p.col << ")" << std::endl;
        }

        pivot p = pivots.back();

        pivots.erase(std::prev(pivots.end()));


        for (unsigned i = 0; i<A.nrow(); ++i) {
            if (i==p.row || A[i].index(p.col)==msparse::row_npos) continue;

            A[i] = row_reduce(p.col, A[i], A[p.row], define_sym);
        }
    }
    std::cout << "done" << std::endl;
}

} // namespace symge
