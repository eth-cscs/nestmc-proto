.. _tutorialsinglecellrecipe:

A simple single cell recipe
===========================

This example builds the same single cell model as
:ref:`the previous tutorial <tutorialsinglecell>`, except using a :class:`arbor.recipe`
and :class:`arbor.simulation` instead of a :class:`arbor.single_cell_model`.

.. Note::

   **Concepts covered in this example:**

   1. Building a :class:`arbor.recipe`.
   2. Using the recipe, default context and domain decomposition to create an :class:`arbor.simulation`
   3. Running the simulation and visualizing the results.

The cell
--------

Let's copy the cell description from the :ref:`previous example <tutorialsinglecell-cell>`,
where construction of the cell is explained in detail.

.. literalinclude:: ../../python/example/single_cell_recipe.py
   :language: python
   :lines: 4,8-23

The recipe
----------

In the :ref:`previous example <tutorialsinglecell-cell>`, the :class:`arbor.single_cell_model` creates
a :class:`arbor.recipe` under the hood, and abstracts away a few details that you may want control over
in more complex simulations. Let's go into those abstractions and create an analogous :class:`arbor.recipe`
manually.

Creating a recipe starts with creating a class that inherits from :class:`arbor.recipe`
and overrides and implements some of :class:`arbor.recipe` methods. Not all methods
have to be overridden, but some will always have to be, such as :meth:`arbor.recipe.num_cells`.
It returns `0` by default and models without cells are quite boring!

.. literalinclude:: ../../python/example/single_cell_recipe.py
   :language: python
   :lines: 25-62

Step **(4)** describes the recipe that will reflect our single cell model.

Step **(4.1)** defines the class constructor. It can take any shape you need, but it
is important to call base class' constructor. If the overridden methods of the class
need to return an object, it may be a good idea to have the returned object be a
member of the class. With this constructor, we could easily change the cell and probes
of the model, should we want to do so. Here we initialize the cell properties to match
Neuron's defaults using Arbor's built-in :meth:`arbor.neuron_cable_properties` and
extend with Arbor's own :meth:`arbor.default_catalogue`.

Step **(4.2)** defines that this model has one cell.

Step **(4.3)** returns :class:`arbor.cell_kind.cable`, the :class:`arbor.cell_kind`
associated with the cable cell defined above. If you mix multiple cell kinds and
descriptions in one recipe, make sure a particular ``gid`` returns matching cell kinds
and descriptions.

Step **(4.4)** returns the cell description passed in on class initialisation. If we
were modelling multiple cells of different kinds, we would need to make sure that the
cell returned by :meth:`arbor.recipe.cell_description` has the same cell kind as
returned by :meth:`arbor.recipe.cell_kind` for every :gen:`gid`.

Step **(4.5)** returns the probes passed in at class initialisation.

Step **(4.6)** returns the properties that will be applied to all cells of that kind in the model.

More methods can be overridden if your model requires that, see :class:`arbor.recipe` for options.

Step **(5)** instantiates the recipe with the cable cell described earlier, and a single voltage probe located at "midpoint".

The context and domain decomposition
------------------------------------

:class:`arbor.single_cell_model` does not only take care of the recipe, it also takes
care of defining how the simulation will be run. When you create and use your own
recipe, you'll need to do this manually, in the form of defining a execution context
and a domain decomposition. Fortunately, the default constructors of
:class:`arbor.context` and :class:`arbor.partition_load_balance` are sufficient for
this model, and is what :class:`arbor.single_cell_model` does under the hood! We'll
leave the details of this subject for another tutorial.

.. literalinclude:: ../../python/example/single_cell_recipe.py
   :language: python
   :lines: 64-67

Step **(6)** sets up a default context and domains.

The simulation
--------------

.. literalinclude:: ../../python/example/single_cell_recipe.py
   :language: python
   :lines: 69-75

Step **(7)** instantiates the simulation and sets up the probe added in step 5. In the
:class:`arbor.single_cell_model` version of this example, the probe frequency and
simulation duration are the same. Note that the frequency is set with a :class:`arbor.regular_schedule`,
which takes a time and not a frequency. Also note that spike recording must be
switched on. For extraction of the probe traces later on, we store a handle.

The results
----------------------------------------------------

Apart from creating :class:`arbor.recipe` ourselves, we have changed nothing
about this simulation compared to :ref:`the previous tutorial <tutorialsinglecell>`.
If we create the same analysis of the results we therefore expect the same results.

.. literalinclude:: ../../python/example/single_cell_recipe.py
   :language: python
   :lines: 77-94

Step **(8)** plots the measured potentials during the runtime of the simulation.
Retrieving the sampled quantities is a little different, these have to be accessed
through the simulation object: :meth:`arbor.simulation.spikes` and :meth:`arbor.simulation.samples`.

We should be seeing something like this:

.. figure:: single_cell_model_result.svg
    :width: 400
    :align: center

    Plot of the potential over time for the voltage probe added in step (5).

You can find the source code for this example in full at ``python/examples/single_cell_recipe.py``.
