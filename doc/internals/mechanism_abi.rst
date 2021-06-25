.. _mechanism_abi:

Mechanism ABI
=============

Here you will find the information needed to connect Arbor to mechanism
implementations outside the use of NMODL and ``modcc``. This may include writing
a custom compiler targetting Arbor, or directly implementing mechanisms in a
C-compatible language. Needless to say that this is aimed at developers rather
than users.

The Arbor library is isolated from these implementations through an Application
Binary Interface (ABI) or plugin interface. Information is provided by the ABI
implementor via two core types.

All functionality is offered via a single C header file in the Arbor include
directory, ``mechanism_abi.h``. The central datatypes here are
``arb_mechanism_type`` and ``arb_mechanism_interface``, laying out the metadata
and backend implementations respectively. A single ``arb_mechanism_type``
instance may be used by multiple ``arb_mechanism_interface`` instances.

Note that ``mechanism_abi.h`` is heavily commented and might be useful as
documentation in its own right.

Metadata: ``arb_mechanism_type``
--------------------------------

This type collects all information independent of the backend.

  .. code:: c

   typedef struct {
     // Metadata
     unsigned long             abi_version;       // mechanism was built using this ABI,
                                                  // should be ARB_MECH_ABI_VERSION
     arb_mechanism_fingerprint fingerprint;       // unique ID, currently ignored
     const char*               name;              // (catalogue-level) unique name
     arb_mechanism_kind        kind;              // one of: point, density, reversal_potential
     bool                      is_linear;         // synapses only: if the state G is governed by dG/dt = f(v, G, M(t)), where:
                                                  //    M(t) =Σ wᵢδᵢ(t) weighted incoming events,
                                                  // then f is linear in G and M. If true, mechanisms must adhere to this contract.
                                                  // Ignored for everything else.
     bool                      has_post_events;   // implements post_event hook
     // Tables
     arb_field_info*           globals;
     arb_size_type             n_globals;
     arb_field_info*           state_vars;
     arb_size_type             n_state_vars;
     arb_field_info*           parameters;
     arb_size_type             n_parameters;
     arb_ion_info*             ions;
     arb_size_type             n_ions;
   } arb_mechanism_type;

Tables
''''''

All tables are given as an integer size and an array. Currently we have two
kinds of tables, which are fairly self-explanatory. Note that these are not
connected to the actual storage layout, in particular, no memory management is
allowed inside mechanisms.

First, parameters, state variables, and global constants

  .. code:: c

    typedef struct {
        const char* name;               // Field name, can be used from the library to query/set field values.
        const char* unit;               // Physical units, just for introspection, not checked
        arb_value_type default_value;   // values will be initialised to this value
        arb_value_type range_low;       // valid range, lower bound, will be enforced
        arb_value_type range_high;      // valid range, upper bound, will be enforced
    } arb_field_info;

Second ion dependencies

  .. code:: c

    typedef struct {
        const char* name;              // Ion name, eg Ca, K, ...
        bool write_int_concentration;  // writes Xi?
        bool write_ext_concentration;  // writes Xo?
        bool write_rev_potential;      // writes Er?
        bool read_rev_potential;       // uses Er?
        bool read_valence;             // Uses valence?
        bool verify_valence;           // Checks valence?
        int  expected_valence;         // Expected value
    } arb_ion_info;

Interlude: Parameter packs
--------------------------

In order to explain the interface type, we have to digress first and introduce
the type ``arb_mechanism_ppack``. This record is used to pass all information to
and from the interface methods.

Objects of this type are always created and allocated by the library and passed
fully formed to the interface. In particular, at this point

- Global data values are initialised
- pointers in ``ion_state_view`` are set to their associated values in shared
  state on the library side
- pointers to state, parameters, globals, and constants are allocated and
  initialised to the given defaults.
- SIMD only: ``index_constraint`` is set up

  .. code:: c

    typedef struct {
        // Global data
        arb_index_type width;                           // Number of CVs of this mechanism, size of arrays
        arb_index_type n_detectors;                     // Number of spike detectors
        arb_index_type* vec_ci;                         // [Array] Map CV to cell
        arb_index_type* vec_di;                         // [Array] Map
        const arb_value_type* vec_t;                    // [Array] time value
        arb_value_type* vec_dt;                         // [Array] time step
        arb_value_type* vec_v;                          // [Array] potential
        arb_value_type* vec_i;                          // [Array] current
        arb_value_type* vec_g;                          // [Array] conductance
        arb_value_type* temperature_degC;               // [Array] Temperature in celsius
        arb_value_type* diam_um;                        // [Array] CV diameter
        arb_value_type* time_since_spike;               // Times since last spike; one entry per cell and detector.
        arb_index_type* node_index;                     // Indices of CVs covered by this mechanism, size is width
        arb_index_type* multiplicity;                   // [Unused]
        arb_value_type* weight;                         // [Array] Weight
        arb_size_type mechanism_id;                     // Unique ID for this mechanism on this cell group
        arb_deliverable_event_stream events;            // Events during the last period
        arb_constraint_partition     index_constraints; // Index restrictions, not initialised for all backends.
        // User data
        arb_value_type** parameters;                    // [Array] setable parameters
        arb_value_type** state_vars;                    // [Array] integrable state
        arb_value_type*  globals;                       // global constant state
        arb_ion_state*   ion_states;                    // [Array] views into shared state
    } arb_mechanism_ppack;

Members tagged as ``[Array]`` represent one value per CV. To access the values
belonging to your mechanism, a level of indirection via ``node_index`` is
needed.

Example: Let's assume mechanism ``hh`` is defined on two regions: ``R``
comprising CVs ``0`` and ``1``, ``R'`` with a single CV ``9``. Then ``node_index
= [0, 1, 9]`` and ``width = 3``. Arrays like ``vec_v`` will be of size ``3`` as
well. To access the CVs' diameters, one would write

  .. code:: c++

    for (auto cv = 0; cv < ppack.width; ++cv) {
        auto idx = node_index[cv];
        auto d   = ppack_um[idx];
    }

Note that values in ``ppack.diam_um`` cover _all_ CV's regardless whether they
are covered by the current mechanisms. Reading those values (or worse writing to
them) is considered undefined behaviour. The same holds for all other fields of
``ppack``.

User Data
'''''''''

This section is derived from the tables passed in via the metadata struct, see
above. One entry per relevant table entry is provided in the same order. So, if

  .. code:: c

    arb_field_info globals[] = { arb_field_info { .name="A",
                                                  .unit="lb ft / s",
                                                  .default_value=42.0,
                                                  .range_low=0,
                                                  .range_high=123 },
                                 arb_field_info { .name="B",
                                                  .unit="kg m / s",
                                                  .default_value=42.0,
                                                  .range_low=0,
                                                  .range_high=123 }};
    arb_mechanism_type m = { .n_globals=2, .globals=globals };

the ``globals`` field of the corresponding parameter pack would have two
entries, the first corresponding to ``A`` and initialised to 42.0 and the second
for ``B`` set to 42.0.

The evolution of the state variables is left to the implementation via
``integrate_state``, while ``globals`` and ``parameters`` are considered
read-only. The ion states internal concentration ``Xi``, external concentration
``Xo``, trans-membrane current ``iX`` may also be read and written. Note that
concurrent updates by multiple mechanisms might occur in any order and each
mechanism will only observe the initial values at the time step boundary. All
contribution by mechanisms are summed up into a final value. Further note that
accessing these values without declaring this via a relevant ``arb_ion_info`` in
the ``arb_mechanism_type`` is undefined behaviour. Parameter packs are specific
to a backend.

Implementation: ``arb_mechanism_interface``
-------------------------------------------

The interface methods will be called with allocated and initialised ``ppack``
data. The actual layout is unspecified, but all pointers are allocated and set
by the library. Plugins should never allocate memory on their own.

  .. code:: C

    typedef void (*arb_mechanism_method)(arb_mechanism_ppack*);

This is the type of all interface methods. These are collected in the record
below with some metadata about the backend.

  .. code:: C

    typedef struct {
      arb_backend_kind   backend;               // one of cpu, gpu
      arb_size_type      partition_width;       // granularity for this backed, eg SIMD lanes
      // Interface methods; see below
      arb_mechanism_method init_mechanism;
      arb_mechanism_method compute_currents;
      arb_mechanism_method apply_events;
      arb_mechanism_method advance_state;
      arb_mechanism_method write_ions;
      arb_mechanism_method post_event;
    } arb_mechanism_interface;


``init_mechanism``
''''''''''''''''''
- called once during instantiation,
- setup initial state, corresponds to NMODL's INITIAL block,
- will receive an allocated and initialised ppack object

``compute_currents``
''''''''''''''''''''

- compute ionic currents and set them through pointers in `ion_state`, currents
  live in `current_density`
- called during each integration time step
  - at the start for reversal potential mechanisms, *before* current reset
  - after event deliver for anything else

``apply_events``
''''''''''''''''

This method is expected to consume a set of `arb_deliverable_events` and apply
effects to internal state, found in ``ppack.events`` which is of type
``arb_deliverable_event_stream``.

  .. code:: c

     typedef struct {
         arb_size_type   mech_id;       // mechanism type identifier (per cell group).
         arb_size_type   mech_index;    // instance of the mechanism
         arb_float_type  weight;        // connection weight
     } arb_deliverable_event;

     typedef struct {
         arb_size_type                n_streams; // number of streams
         const arb_deliverable_event* events;    // array of event data items
         const arb_index_type*        begin;     // array of offsets to beginning of marked events
         const arb_index_type*        end;       // array of offsets to end of marked events
     }  arb_deliverable_event_stream;

These structures are set up correctly externally, but are only valid during this call.
The data is read-only for ``apply_events``.

- called during each integration time step, right after resetting currents
- corresponding to ``NET_RECEIVE``

``advanced_state``
''''''''''''''''''

- called during each integration time step, after solving Hines matrices
- perform integration on state variables
- state variables live in `state_vars`, with a layout described above

``write_ions``
''''''''''''''

- update ionic concentrations via the pointers in `ion_state`
- called during each integration time step, after state integration

``post_event``
''''''''''''''

- used to implement spike time dependent plasticity
- consumes ``ppack.time_since_spike``
- called during each integration time step, after checking for spikes
- if implementing this, also set ``has_post_event=true`` in the metadata

SIMDization
-----------

If a mechanism interface processes arrays in SIMD bundles, it needs to set
``partition_width`` to that bundle's width in units of ``arb_value_type``. The
library will set up ``arb_constraint_partition index_constraint`` in the
parameter pack. This structure describe which bundles can be loaded/stored as a
contiguous block, which ones must be gathered/scattered, which are to be
broadcast from a constant, and so on. The reason for this is the indirection via
``node_index`` mentioned before. Please refer to the documentation of our SIMD
interface layer for more information.

Making A Loadable Mechanism
---------------------------

Mechanisms interface with the library by providing three functions, one
returning the metadata portion, and one for each implemented backend (currently
two). The latter may return a NULL pointer, indicating that this backend is not
supported. The naming scheme is shown in the example below

  .. code:: C

    arb_mechanism_type make_arb_default_catalogue_pas();

    arb_mechanism_interface* make_arb_default_catalogue_pas_interface_multicore();
    arb_mechanism_interface* make_arb_default_catalogue_pas_interface_gpu();
