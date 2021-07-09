.. _tutorialsinglecell:

A simple single cell model
==========================

Building and testing detailed models of individual cells, then optimizing their
parameters is usually the first step in building models with multi-compartment cells.
Arbor supports a *single cell model* workflow for this purpose, which is a good way to
introduce Arbor's cell modelling concepts and approach.

.. Note::

   **Concepts covered in this example:**

   1. Intro to building a morphology from a :class:`arbor.segment_tree`.
   2. Intro to region and locset expressions.
   3. Intro to decors and cell decorations.
   4. Building a :class:`arbor.cable_cell` object.
   5. Building a :class:`arbor.single_cell_model` object.
   6. Running a simulation and visualising the results.

.. _tutorialsinglecell-cell:

The cell
--------

The most trivial representation of a cell in Arbor is to model the entire cell as a
single cylinder. The following example shows the steps required to construct a model of a
cylindrical cell with a length of 6 μm and a radius of 3 μm; Hodgkin–Huxley dynamics
and a current clamp stimulus, then run the model for 30 ms.

The first step is to construct the cell. In Arbor, the abstract representation used to
define a cell with branching cable morphology is a ``cable_cell``, which holds a
description of the cell's morphology, named regions and locations on the morphology, and
descriptions of ion channels, synapses, spike detectors and electrical properties.

Our *single-segment HH cell* has a simple morphology and dynamics, constructed as follows:

.. literalinclude:: ../../python/example/single_cell_model.py
   :language: python
   :lines: 4,6-23

Step **(1)** constructs a :class:`arbor.segment_tree` (see also :ref:`segment tree<morph-segment_tree>`).
The segment tree is the representation used to construct the morphology of a cell. A segment is
a tapered cone with a tag; the tag can be used to classify the type of the segment (for example
soma, dendrite etc). To create a segment tree representing our single-cylinder cell, we need to add
one segment to our ``tree`` object. We use the :meth:`arbor.segment_tree.append` method, which takes
4 arguments: the parent segment which does not exist for the first segment, so we use :class:`arbor.mnpos`;
the proximal :class:`arbor.mpoint` (location and radius) of the segment; the distal :class:`arbor.mpoint`
of the segment; and the tag.

Step **(2)** creates a dictionary of labels (:class:`arbor.label_dict<arbor.label_dict>`). Labels give
names to :term:`regions<region>` and :term:`location<locset>` described using a DSL
based on s-expressions. Labels from the dictionary can then be used to facilitate adding synapses,
dynamics, stimuli and probes to the cell. We add two labels:

* ``soma`` defines a *region* with ``(tag  1)``. Note that this corresponds to the
  ``tag`` parameter that was used to define the single segment in step (1).
* ``midpoint`` defines a *location* at ``(location 0 0.5)``, which is the mid point ``0.5``
  of branch ``0``, which corresponds to the midpoint of the soma on the morphology defined in step (1).

Step **(3)** constructs a :class:`arbor.decor` that describes the distribution and placement
of dynamics and properties on a cell.  The cell's default properties can be modified, and we can use
:meth:`arbor.decor.paint` and :meth:`arbor.decor.place` to further customise it in the
following way:

* :meth:`arbor.decor.set_property` is used to set some default properties on the entire cell.
  In the above example we set the initial membrane potential to -40 mV.
* :meth:`arbor.decor.paint` is used to set properties or add dynamics to a region of the cell.
  We call this method 'painting' to convey that we are working on sections of a cell, as opposed to
  precise locations: for example, we might want to *paint* an ion channel on all dendrites, and then
  *place* a synapse at the tip of the axon. In the above example we paint
  HH dynamics on the region we previously named ``"soma"`` in our label dictionary.
* :meth:`arbor.decor.place` is used to add objects on a precise
  :class:`arbor.location` on a cell. Examples of objects that are *placed* are synapses,
  spike detectors, current stimuli, and probes. In the above example we place a current stimulus
  :class:`arbor.iclamp` with a duration of 2 ms and a current of 0.8 nA, starting at 10 ms
  on the location we previously labelled ``"midpoint"``. We also place a :class:`arbor.spike_detector`
  with a threshold of -10 mV on the same location.

Step **(4)** constructs the :class:`arbor.cable_cell` from the segment tree and dictionary of labelled regions and locations.

The single cell model
---------------------

Once the cell description has been built, the next step is to build and run the simulation.
Arbor provides an interface for constructing single cell models with the
:class:`arbor.single_cell_model` helper that creates a model from a cell description, with
an interface for recording outputs and running the simulation.

The single cell model has 4 main functions:

1. It holds the **global properties** of the model
2. It registers **probes** on specific locations on the cell to measure the voltage.
3. It **runs** the simulation.
4. It collects **spikes** from spike detectors and voltage **traces** from registered probes.

Right now, we'll only set a probe and run the simulation.

.. literalinclude:: ../../python/example/single_cell_model.py
   :language: python
   :lines: 25-32

Step **(5)** instantiates the :class:`arbor.single_cell_model`
with our single-compartment cell.

Step **(6)** adds a :meth:`arbor.single_cell_model.probe`
used to record variables from the model. Three pieces of information are
provided: the type of quantity we want probed (voltage), the location where we want to
probe ('"midpoint"'), and the frequency at which we want to sample (10 kHz).

Step **(7)** runs the actual simulation for a duration of 30 ms.

The results
-----------

Our cell and model have been defined and we have run our simulation. Now we can look at what
the spike detector and a voltage probes from our model have produced.

.. literalinclude:: ../../python/example/single_cell_model.py
   :language: python
   :lines: 34-46

Step **(8)** accesses :meth:`arbor.single_cell_model.spikes`
to print the spike times. A single spike should be generated at around the same time the stimulus
we provided in step (3) gets activated (10ms).

Step **(9)** plots the measured potentials during the runtime of the simulation. The sampled quantities
can be accessed through :meth:`arbor.single_cell_model.traces`.
We should be seeing something like this:

.. figure:: single_cell_model_result.svg
    :width: 400
    :align: center

    Plot of the potential over time for the voltage probe added in step (6).

The full code
-------------

You can find the source code for this example in full at ``python/examples/single_cell_model.py``.
