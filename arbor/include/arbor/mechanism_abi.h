#ifndef ARB_MECH_ABI
#define ARB_MECH_ABI

#include <arbor/arb_types.h>

#ifdef __cplusplus
extern "C" {
#endif

// Version
#define ARB_MECH_ABI_VERSION_MAJOR 0
#define ARB_MECH_ABI_VERSION_MINOR 0
#define ARB_MECH_ABI_VERSION_PATCH 1
#define ARB_MECH_ABI_VERSION ((ARB_MECH_ABI_VERSION_MAJOR * 10000 * 10000) + (ARB_MECH_ABI_VERSION_MAJOR * 10000) + ARB_MECH_ABI_VERSION_PATCH)

typedef const char* arb_mechanism_fingerprint;

// Selectors
typedef uint32_t arb_mechanism_kind;
#define arb_mechanism_kind_nil 0
#define arb_mechanism_kind_point 1
#define arb_mechanism_kind_density 2
#define arb_mechanism_kind_reversal_potential 3

typedef uint32_t arb_backend_kind;
#define arb_backend_kind_nil 0
#define arb_backend_kind_cpu 1
#define arb_backend_kind_gpu 2

// Ion state variables; view into shared_state
typedef struct arb_ion_state {
    arb_value_type* current_density;
    arb_value_type* reversal_potential;
    arb_value_type* internal_concentration;
    arb_value_type* external_concentration;
    arb_value_type* ionic_charge;
    arb_index_type* index;
} arb_ion_state;

// Event; consumed by `apply_event`
typedef struct arb_deliverable_event_data {
    arb_size_type   mech_id;       // mechanism type identifier (per cell group).
    arb_size_type   mech_index;    // instance of the mechanism
    arb_weight_type weight;
} arb_deliverable_event_data;

/* A set of `n` streams of events, where those in the
 * ranges (events + begin[i], events + end[i]) i = 0..n-1
 * are meant to be consumed
 */
typedef struct arb_deliverable_event_stream {
    arb_size_type                     n_streams; // number of streams
    const arb_deliverable_event_data* events;    // array of event data items
    const arb_index_type*             begin;     // array of offsets to beginning of marked events
    const arb_index_type*             end;       // array of offsets to end of marked events
}  arb_deliverable_event_stream;

// Constraints for use in SIMD implementations, see there.
typedef struct arb_constraint_partition {
    arb_size_type   n_contiguous;
    arb_size_type   n_constant;
    arb_size_type   n_independent;
    arb_size_type   n_none;
    arb_index_type* contiguous;
    arb_index_type* constant;
    arb_index_type* independent;
    arb_index_type* none;
} arb_constraint_partition;

// Parameter Pack
typedef struct arb_mechanism_ppack {
    arb_size_type  width;                           // Number of CVs
    arb_index_type  n_detectors;                     // Number of spike detectors
    arb_index_type* vec_ci;
    arb_index_type* vec_di;
    const arb_value_type* vec_t;
    arb_value_type* vec_dt;
    arb_value_type* vec_v;
    arb_value_type* vec_i;
    arb_value_type* vec_g;
    arb_value_type* temperature_degC;
    arb_value_type* diam_um;
    arb_value_type* time_since_spike;
    arb_index_type* node_index;
    arb_index_type* multiplicity;
    arb_value_type* weight;
    arb_size_type   mechanism_id;

    arb_deliverable_event_stream events;            // Events during the last period
    arb_constraint_partition     index_constraints; // Index restrictions, not initialised for all backend.

    arb_value_type** parameters;                    // Array of setable parameters     (Array)
    arb_value_type** state_vars;                    // Array of integrable state       (Array)
    arb_value_type*  globals;                       // Array of global constant state  (Scalar)
    arb_ion_state*   ion_states;                    // Array of views into shared state
} arb_mechanism_ppack;


/* Mechanism Plugin
 *
 * Everything below has to be filled out by the plugin author/compiler.
 * The interface methods will be called with allocated and initialised `ppack`
 * data. The actual layout is unspecified, but all pointers are allocated and set
 * by the library. Plugins should never allocate memory on their own.
 */
typedef void (*arb_mechanism_method)(arb_mechanism_ppack*); // Convenience for extension methods

typedef struct arb_mechanism_interface {
    arb_backend_kind   backend;               // GPU, CPU, ...
    arb_size_type      partition_width;       // Width for partitioning indices, based on SIMD for example
    arb_size_type      alignment;
    // Interface methods; hooks called by the engine during the lifetime of the mechanism.
    /* 1. init_mechanism
     *   - called once during instantiation,
     *   - setup initial state, corresponds to NMODL's INITIAL block,
     *   - will receive an allocated and initialised ppack object
     *     - pointers in ion_state_view are set to their associated values in shared state
     *     - pointers to state, parameters, globals, and constants are allocated and initialised to the given defaults.
     *     - SIMD only: index_constraint is set up
     *     - Internal values (see above) are initialised
     */
    arb_mechanism_method init_mechanism;
    /* 2. compute_currents
     * - compute ionic currents
     *   - pointers in `ion_state` are set to [ion_0, ion_1, ...] from the `ions` table
     *   - currents live in `current_density`
     * - called during each integration time step
     *   - at the start for reversal potential mechanisms, *before* current reset
     *   - after event deliver for anything else
     */
    arb_mechanism_method compute_currents;

    /* 3. apply_events
     * - consume `deliverable_events` and apply effects to internal state
     * - `deliverable_events` is setup correctly externally, is read-only for apply events
     * - called during each integration time step, right after resetting currents
     */
    arb_mechanism_method apply_events;
    /* 4. advanced_state
     * - called during each integration time step, after solving Hines matrices
     * - perform integration on state variables, often given as an ODE
     * - state variables live in `ppack::state_vars`
     */
    arb_mechanism_method advance_state;
    /* 5. write_ions
     * - update ionic concentrations
     *   - pointers in `ion_state` are set to [ion_0, ion_1, ...] from the `ions` table
     *   - variables live in `internal_concentration` and `external_concentration`
     * - called during each integration time step, after state integration
     */
    arb_mechanism_method write_ions;
    /* 6. post_event
     * - called during each integration time step, after checking for spikes
     * - corresponds to NET_RECEIVE in NMODL
     */
    arb_mechanism_method post_event;
} arb_mechanism_interface;

typedef struct arb_field_info {
    const char* name;
    const char* unit;
    arb_value_type default_value;
    arb_value_type range_low;
    arb_value_type range_high;
} arb_field_info;

// Ion dependency
typedef struct arb_ion_info {
    const char* name;
    bool write_int_concentration;
    bool write_ext_concentration;
    bool write_rev_potential;
    bool read_rev_potential;
    bool read_valence;
    bool verify_valence;
    int  expected_valence;
} arb_ion_info;

// Backend independent data
typedef struct arb_mechanism_type {
    // Metadata
    unsigned long             abi_version;      // plugin ABI version used to build this mechanism
    arb_mechanism_fingerprint fingerprint;      // provide a unique ID
    const char*               name;             // provide unique name
    arb_mechanism_kind        kind;             // Point, Density, ReversalPotential, ...
    bool                      is_linear;        // linear, homogeneous mechanism
    bool                      has_post_events;
    // Tables
    arb_field_info*           globals;          // Global constants
    arb_size_type             n_globals;
    arb_field_info*           state_vars;       // Integrable state
    arb_size_type             n_state_vars;
    arb_field_info*           parameters;       // Mechanism parameters
    arb_size_type             n_parameters;
    arb_ion_info*             ions;             // Ion properties
    arb_size_type             n_ions;
} arb_mechanism_type;

#ifdef __cplusplus
}
#endif
#endif
