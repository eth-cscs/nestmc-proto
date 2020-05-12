#pragma once

#include <cstddef>
#include <memory>
#include <unordered_map>
#include <stdexcept>

#include <arbor/arbexcept.hpp>
#include <arbor/common_types.hpp>
#include <arbor/event_generator.hpp>
#include <arbor/util/unique_any.hpp>

namespace arb {

struct probe_info {
    cell_member_type id;
    probe_tag tag;

    // Address type will be specific to cell kind of cell `id.gid`.
    util::any address;
};

/* Recipe descriptions are cell-oriented: in order that the building
 * phase can be distributed, and in order that the recipe description
 * can be built indepedently of any runtime execution environment.
 */

// Note: `cell_connection` and `connection` have essentially the same data
// and represent the same thing conceptually. `cell_connection` objects
// are notionally described in terms of external cell identifiers instead
// of internal gids, but we are not making the distinction between the
// two in the current code. These two types could well be merged.

struct cell_connection {
    // Connection end-points are represented by pairs
    // (cell index, source/target index on cell).
    using cell_connection_endpoint = cell_member_type;

    cell_connection_endpoint source;
    cell_connection_endpoint dest;

    float weight;
    float delay;

    cell_connection(cell_connection_endpoint src, cell_connection_endpoint dst, float w, float d):
        source(src), dest(dst), weight(w), delay(d)
    {}
};

struct gap_junction_connection {
    cell_member_type local;
    cell_member_type peer;
    double ggap;

    gap_junction_connection(cell_member_type local, cell_member_type peer, double g):
            local(local), peer(peer), ggap(g) {}
};

class recipe {
public:
    virtual cell_size_type num_cells() const = 0;

    // Cell description type will be specific to cell kind of cell with given gid.
    virtual util::unique_any get_cell_description(cell_gid_type gid) const = 0;
    virtual cell_kind get_cell_kind(cell_gid_type) const = 0;

    virtual cell_size_type num_sources(cell_gid_type) const { return 0; }
    virtual cell_size_type num_targets(cell_gid_type) const { return 0; }
    virtual cell_size_type num_probes(cell_gid_type)  const { return 0; }
    virtual cell_size_type num_gap_junction_sites(cell_gid_type gid)  const {
        return gap_junctions_on(gid).size();
    }
    virtual std::vector<event_generator> event_generators(cell_gid_type) const {
        return {};
    }
    virtual std::vector<cell_connection> connections_on(cell_gid_type) const {
        return {};
    }
    virtual std::vector<gap_junction_connection> gap_junctions_on(cell_gid_type) const {
        return {};
    }

    virtual probe_info get_probe(cell_member_type probe_id) const {
        throw bad_probe_id(probe_id);
    }

    // Global property type will be specific to given cell kind.
    virtual util::any get_global_properties(cell_kind) const { return util::any{}; };

    virtual ~recipe() {}
};

} // namespace arb

static inline std::ostream& operator<<(std::ostream& o, arb::cell_connection const& con) {
    return o << "ccon [" << con.source << " -> " << con.dest
             << " : weight " << con.weight
             << ", delay " << con.delay << "]";
}
