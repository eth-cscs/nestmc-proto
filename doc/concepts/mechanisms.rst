.. _mechanisms:

Cable cell mechanisms
=====================

Mechanisms describe biophysical processes such as ion channels and synapses.
Mechanisms are assigned to regions and locations on a cell morphology
through the process of :ref:`decoration <cablecell-decoration>`.
Mechanisms are described using a dialect of the :ref:`NMODL <nmodl>` domain
specific language that is similarly used in `NEURON <https://neuron.yale.edu/neuron/>`_.

Arbor supports mechanism descriptions using the NMODL language through our ``modcc``
compiler. ``modcc`` supports many of NMODL's features but there are a few
additional :ref:`guidelines <formatnmodl>`.
for users who wish to compile their own mechanisms for Arbor. Unfortunately, out-of-tree
mechanism building is not yet supported in Arbor. However, we have many built-in mechanisms
that can be used, which are oragnized in *mechanism catalogues*.

Mechanism catalogues
--------------------

A *mechanism catalogue* is a collection of mechanisms that maintains:

1. A collection of mechanism metadata indexed by name.
2. A further hierarchy of *derived* mechanisms, that allow specialization of
   global parameters, ion bindings, and implementations.
3. A map for looking up a concrete mechanism implementation on a target hardware back end.

Derived mechanisms will always have a different name to the mechanism from which they are derived.
This name is given explicitly when the derivation is constructed, or implicitly when a mechanism
is :ref:`requested <mechanisms-name-note>` with a name of the form ``"mech/param=value,..."``.
In this instance, if a mechanism of that name does not already exist in the catalogue, it will be
implicitly derived from an existing mechanism ``"mech"`` with global parameters and ion bindings
set according to the assignments following the slash. If the mechanism ``"mech"`` depends upon
only a single ion, the name of that ion can be omitted in the assignments:
``"mech/oldion=newion"`` and ``"mech/newion"`` are equivalent derivations.


Catalogues provide an interface for querying mechanism metadata, which includes the following information:

* Global parameter names, units, and default values.
* Range parameter names, units, and default values.
* State variable names, units and default values.
* Ion dependencies: for each ion used by the mechanism, information on whether the mechanism writes
  to its internal or external concentration or to its reversal potential value, and whether it reads
  or asserts the ionic charge.

Arbor provides a default catalogue of mechanisms as well as two other catalogues containing the sets of common mechanisms
used by the `Allen Institute <https://alleninstitute.org/>`_ and the `Blue Brain Project <https://portal.bluebrain.epfl.ch/>`_.
(Find the NMODL descriptions of the `default mechanisms <https://github.com/arbor-sim/arbor/tree/master/mechanisms/default>`_,
the `Allen institute mechanisms <https://github.com/arbor-sim/arbor/tree/master/mechanisms/allen>`_ and
the `BBP mechanisms <https://github.com/arbor-sim/arbor/tree/master/mechanisms/bbp>`_ at the provided links.)

.. _mechanisms_builtins:

Built-in Catalogues
'''''''''''''''''''

Arbor provides the *default* catalogue with the following mechanisms:

* *pas*: Leaky current (:ref:`density mechanism <mechanisms-density>`).
* *hh*: Classic Hodgkin-Huxley dynamics (:ref:`density mechanism
  <mechanisms-density>`).
* *nernst*: Calculate reversal potential for an ionic species using the Nernst
  equation (:ref:`reversal potential mechanism <mechanisms-revpot>`). **NB**
  This is not meant to be used directly
* *expsyn*: Synapse with discontinuous change in conductance at an event
  followed by an exponential decay (:ref:`point mechanism <mechanisms-point>`).
* *exp2syn*: Bi-exponential conductance synapse described by two time constants:
  rise and decay (:ref:`point mechanism <mechanisms-point>`).

With the exception of *nernst*, these mechanisms are the same as those available in NEURON.

Two catalogues are provided that collect mechanisms associated with specific projects and model databases:

* *bbp* For models published by the Blue Brain Project (BBP).
* *allen* For models published on the Allen Brain Atlas Database.

.. _mechanisms_dynamic:

Adding Catalogues to Arbor
''''''''''''''''''''''''''

.. Note::

   If you are coming from NEURON this is the equivalent of ``nrnivmodl``.

This will produce a catalogue loadable at runtime by calling ``load_catalogue``
with a filename in both C++ and Python. The steps are

1. Prepare a directory containing your NMODL files (.mod suffixes required)
2. Call ``build-catalogue`` installed by arbor

   .. code-block :: bash

     build-catalogue <name> <path/to/nmodl>

All files with the suffix ``.mod`` located in ``<path/to/nmodl>`` will be baked into
a catalogue named ``lib<name>-catalogue.so`` and placed into your current working
directory. Note that these files are platform-specific and should only be used
on the combination of OS, compiler, arbor, and machine they were built with.

Errors might be diagnosable by passing the ``-v`` flag.

This catalogue can then be used similarly to the built-in ones

   .. code-block :: python

     import arbor as A

     c = A.load_catalogue('bbp2-catalogue.so')

     [n for n in c]
     >> ['Ca_LVAst',
         'Nap_Et2',
         'NaTa_t',
         'SKv3_1',
         'K_Tst',
         'Ih',
         'SK_E2',
         'Ca_HVA',
         'CaDynamics_E2',
         'Im',
         'NaTs2_t',
         'K_Pst']

See also the demonstration in ``python/example/dynamic-catalogue.py`` for an example.

Parameters
''''''''''

Mechanism behaviour can be tuned using parameters and ion channel dependencies,
as defined in the NMODL description.
Parameters and ion species are set initially before a simulation starts, and remain
unchanged thereafter, for the duration of the simulation.
There are two types of parameters that can be set by users:

* *Global* parameters are a single scalar value that is the same everywhere a mechanism is defined.
* *Range* parameters can vary spatially.

Every mechanism is applied to a cell via a *mechanism description*, a
``(name, range_parameters)`` tuple, where ``name`` is a string,
and ``range_parameters`` is an optional dictionary of key-value pairs
that specifies values for range parameters.
For example, consider a mechanism that models passive leaky dynamics with
the following parameters:

* *Name*: ``"pas"``.
* *Global parameter*: reversal potential ``e``, default -70 mV.
* *Range parameter*: conductance ``g``, default 0.001 S⋅cm⁻².

The following example mechanism descriptions for our passive mechanism show that parameters and
ion species dependencies only need to be specified when they differ from their defaults:

* ``("pas")``: the passive mechanism with default parameters.
* ``("pas/e=-80")``: derive a new passive mechanism with a non-default value for global parameter.
* ``("pas", {"g": 0.005})``: passive mechanism with a new a non-default range parameter value.
* ``("pas/e=-80", {"g": 0.005})``: derive a new passive mechanism that overrides both

Similarly to global parameters, ion species can be renamed in the mechanism name.
This allows the use of generic mechanisms that can be adapted to a specific species
during model instantiation.
For example, the ``nernst`` mechanism in Arbor's default mechanism catalogue calculates
the reversal potential of a generic ionic species ``x`` according to its internal
and external concentrations and valence. To specialize ``nernst`` for calcium name it
``("nernst/x=ca")``, or as there is only one ion species in the mechanism the
shorthand ``("nernst/ca")`` can be used unambiguously.

.. _mechanisms-name-note:

.. note::
    Global parameter values and ionic dependencies are the same for each instance of
    a mechanism; changing these requires the derivation of a new mechanism, implicitly or explicitly.
    For this reason, new global parameter values and ion renaming are part of the name of
    the new mechanism, or a mechanism with a new unique name must be defined.


Mechanism types
---------------

There are two broad categories of mechanism, density mechanisms and
point mechanisms, and a third special density mechanism for
computing ionic reversal potentials.

.. _mechanisms-density:

Density mechanisms
''''''''''''''''''''''

Density mechanisms describe biophysical processes that are distributed in space,
but whose behaviour is defined purely by the state of the cell and the process
at any given point.

Density mechanisms are commonly used to describe ion channel dynamics,
for example the ``hh`` and ``pas`` mechanisms provided by NEURON and Arbor,
which model classic Hodgkin-Huxley and passive leaky currents respectively.

.. _mechanisms-revpot:

Ion reversal potential mechanisms
'''''''''''''''''''''''''''''''''

These mechanisms, which describe ionic reversal potential
behaviour, can be specified for cells or the whole model.

The reversal potential of an ion species is calculated by an
optional *reversal potential mechanism*.
If no such mechanism is specified for an ion species, the initial
reversal potential values are maintained for the course of a simulation.
Otherwise, the mechanism does the work.

Reversal potential mechanisms are density mechanisms subject to some strict restrictions.
Specifically, a reversal potential mechanism described in NMODL:

* May not maintain any state variables.
* Can only write to the reversal potential (``eX``) value of the ion species.
* Can not be a :ref:`point mechanism <mechanisms-point>`.

Essentially, reversal potential mechanisms must be pure functions of cellular
and ionic state.

.. note::
    Arbor imposes greater restrictions on mechanisms that update ionic reversal potentials
    than NEURON. Doing so simplifies reasoning about interactions between
    mechanisms that share ionic species, by virtue of having one mechanism, and one
    mechanism only, that calculates reversal potentials according to concentrations
    that the other mechanisms use and modify.

.. _mechanisms-point:

Point mechanisms
'''''''''''''''''''''''''''''''''

*Point mechanisms*, which are associated with connection end points on a
cable cell, are placed at discrete locations on the cell.
Unlike density mechanisms, whose behaviour is defined purely by the state of the cell
and the process, their behaviour is additionally governed by the timing and weight of
events delivered via incoming connections.


API
---

* :ref:`Python <py_mechanisms>`
* *TODO* C++ documentation.
