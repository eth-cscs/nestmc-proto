.. _tutorialsinglecellswcrecipe:

A detailed single cell recipe
=============================

This example builds the same single cell model as
:ref:`the previous tutorial <tutorialsinglecellswc>`, except using a :class:`arbor.recipe`
and :class:`arbor.simulation` instead of a :class:`arbor.single_cell_model`.

.. Note::

   **Concepts covered in this example:**

   1. Building a :class:`arbor.recipe`.
   2. Building an :class:`arbor.context` and a :class:`arbor.domain_decomposition`
   3. Using the recipe, context and domain decomposition to create a :class:`arbor.simulation`
   4. Running the simulation and visualizing the results,


Recipes are an important concept in Arbor. They represent the most versatile tool
for building a complex network of cells. We will go though this example of a model
of a single cell, before using the recipe to represent more complex networks in
subsequent examples.

We outline the following steps of this example:

1. Define the **cell**. This is the same cell we have seen before.
2. Define the **recipe** of the model.
3. Define the **execution context** of the model: a description of the underlying system
   on which the simulation will run.
4. Define the **domain decomposition** of the network: how the cells are distributed on
   the different ranks of the system.
5. Define the **simulation**.
6. **Run** the simulation.
7. Collect and visualize the **results**.

The cell
********

We can copy the cell description code or reuse ``single_cell_detailed.swc`` from the
:ref:`previous example <tutorialsinglecellswc-cell>` where it is explained in detail.

We will need to add one more thing to the cell. We will create the voltage probe at the "custom_terminal" locset.
In the previous example, this probe was registered directly using the :class:`arbor.single_cell_model` object.
Now it has to be explicitly created and registered in the recipe.

.. _tutorialsinglecellswcrecipe-probe:

.. literalinclude:: ../../python/example/single_cell_detailed_recipe.py
   :language: python
   :lines: 87-89

The recipe
**********

The :class:`arbor.single_cell_model` of the previous example created a :class:`arbor.recipe` under
the hood, and abstracted away the details so we were unaware of its existence. In this example, we will
examine the recipe in detail: how to create one, and why it is needed.

.. literalinclude:: ../../python/example/single_cell_detailed_recipe.py
   :language: python
   :lines: 91-143

Let's go through the recipe point by point.

Step **(6)** creates a ``single_recipe`` class that inherits from :class:`arbor.recipe`. The base recipe
implements all the methods defined above with default values except :meth:`arbor.recipe.num_cells`,
:meth:`arbor.recipe.cell_kind` and :meth:`arbor.recipe.cell_description` which always have to be implemented
by the user. The :meth:`arbor.recipe.global_properties` also needs to be implemented for
:class:`arbor.cell_kind.cable` cells. The inherited recipe can implement any number of additional methods and
have any number of instance or class variables.

Step **(6.1)** defines the class constructor. In this case, we pass a ``cell`` and a set of ``probes`` as
arguments. These will be used to initialize the instance variables ``self.the_cell`` and ``self.the_probes``,
which will be used in the overloaded ``cell_description`` and ``get_probes`` methods. Before variable
initialization, we call the base C++ class constructor ``arbor.recipe.__init__(self)``. This ensures correct
initialization of memory in the C++ class.

We also create the ``self.the_cat`` variable and set it to arbor's default mechanism catalogue. This will expose
the *hh* and *pas* mechanisms but not the *Ih* mechanism, which is present in the allen catalogue. To be able
to use *Ih*, we extend ``self.the_cat`` to include the allen catalogue.

Finally we create the ``self.the_props`` variable. This will hold the global properties of the model, which apply
to all the cells in the network. Initially it is empty. We set all the properties of the system similar to
what we did in the :ref:`previous example <tutorialsinglecellswc-gprop>`. One last important step is to register
``self.the_cat`` with ``self.the_props``.

.. Note::

   The mechanism catalogue needs to live in the recipe as an instance variable. Its lifetime needs to extend
   to the entire duration of the simulation.

Step **(6.2)** overrides the :meth:`arbor.recipe.num_cells` method. It takes no arguments. We simply return 1,
as we are only simulating one cell in this example.

Step **(6.3)** overrides the :meth:`arbor.recipe.cell_kind` method. It takes one argument: ``gid``.
Given the gid, this method returns the kind of the cell. Our defined cell is a
:class:`arbor.cell_kind.cable`, so we simply return that.

Step **(6.4)** overrides the :meth:`arbor.recipe.cell_description` method. It takes one argument: ``gid``.
Given the gid, this method returns the cell description which is the cell object passed to the constructor
of the recipe. We return ``self.the_cell``.

Step **(6.5)** overrides the :meth:`arbor.recipe.get_probes` method. It takes one argument: ``gid``.
Given the gid, this method returns all the probes on the cell. The probes can be of many different kinds
measuring different quantities on different locations of the cell. We pass these probes explicitly to the recipe
and they are stored in ``self.the_probes``, so we return that variable.

Step **(6.6)** overrides the :meth:`arbor.recipe.connections_on` method. It takes one argument: ``gid``.
Given the gid, this method returns all the connections ending on that cell. These are typically synapse
connections from other cell *sources* to specific *targets* on the cell with id ``gid``. Since we are
simulating a single cell, and self-connections are not possible, we return an empty list.

Step **(6.7)** overrides the :meth:`arbor.recipe.gap_junctions_on` method. It takes one argument: ``gid``.
Given the gid, this method returns all the gap junctions on that cell. Gap junctions require 2 separate cells.
Since we are simulating a single cell, we return an empty list.

Step **(6.8)** overrides the :meth:`arbor.recipe.event_generators` method. It takes one argument: ``gid``.
Given the gid, this method returns *event generators* on that cell. These generators trigger events (or
spikes) on specific *targets* on the cell. They can be used to simulate spikes from other cells, to kick-start
a simulation for example. Our cell uses a current clamp as a stimulus, and has no targets, so we return
an empty list.

Step **(6.9)** overrides the :meth:`arbor.recipe.global_properties` method. It takes one argument: ``kind``.
This method returns the default global properties of the model which apply to all cells in the network of
that kind. We return ``self.the_props`` which we defined in step **(1)**.

.. Note::

   You may wonder why the method :meth:`arbor.recipe.cell_kind` is required, since it can be inferred by examining the cell description.
   The recipe was designed to allow building simulations efficiently in a distributed system with minimum
   communication. Some parts of the model initialization require only the cell kind,
   not the full cell description which can be quite expensive to build. Providing these
   descriptions separately saves time and resources for the user.

   More information on the recipe can be found :ref:`here <modelrecipe>`.

Now we can instantiate a ``single_recipe`` object using the ``cell`` and ``probe`` we created in the
previous section:

.. literalinclude:: ../../python/example/single_cell_detailed_recipe.py
   :language: python
   :lines: 145-147

The execution context
*********************

An :ref:`execution context <modelcontext>` describes the hardware resources on which the simulation will run.
It contains the thread pool used to parallelise work on the local CPU, and optionally describes GPU resources
and the MPI communicator for distributed simulations. In the previous
examples, the :class:`arbor.single_cell_model` object created the execution context :class:`arbor.context`
behind the scenes.

The details of the execution context can be customized by the user. We may specify the number of threads
in the thread pool; determine the id of the GPU to be used; or create our own MPI communicator. However,
the ideal settings can usually be inferred from the system, and Arbor can do that with a simple command.

.. literalinclude:: ../../python/example/single_cell_detailed_recipe.py
   :language: python
   :lines: 149-151

The domain decomposition
************************

The domain decomposition describes the distribution of the cells over the available computational resources.
The :class:`arbor.single_cell_model` also handled that without our knowledge in the previous examples.
Now, we have to define it ourselves.

The :class:`arbor.domain_decomposition` class can be manually created by the user, by deciding which cells
go on which ranks. Or we can use a load balancer that can partition the cells across ranks according to
some rules. Arbor provides :class:`arbor.partition_load_balance`, which, using the recipe and execution
context, creates the :class:`arbor.domain_decomposition` object for us.

Our example is a simple one, with just one cell. We don't need any sophisticated partitioning algorithms, so
we can use the load balancer, which does a good job distributing simple networks.

.. literalinclude:: ../../python/example/single_cell_detailed_recipe.py
   :language: python
   :lines: 153-155

The simulation
**************

Finally we have the 3 components needed to create a :class:`arbor.simulation` object.

.. literalinclude:: ../../python/example/single_cell_detailed_recipe.py
   :language: python
   :lines: 157-159

Before we run the simulation, however, we need to register what results we expect once execution is over.
This was handled by the :class:`arbor.single_cell_model` object in the previous example.

We would like to get a list of the spikes on the cell during the runtime of the simulation, and we would like
to plot the voltage registered by the probe on the "custom_terminal" locset.

.. literalinclude:: ../../python/example/single_cell_detailed_recipe.py
   :language: python
   :lines: 161-166

The lines handling probe sampling warrant a second look. First, we declared ``probe_id`` to be a
:class:`arbor.cell_member`, with :class:`arbor.cell_member.gid` = 0 and :class:`arbor.cell_member.index` = 0.
This variable serves as a global identifier of a probe on a cell, namely the first declared probe on the
cell with gid = 0, which is id of the :ref:`only probe <tutorialsinglecellswcrecipe-probe>` we created on
the only cell in the model.

Next, we instructed the simulation to sample ``probe_id`` at a frequency of 50 kHz. That function returns a
``handle`` which we will use to extract the results of the sampling after running the simulation.

The execution
*************

We can now run the simulation we just instantiated for a duration of 100 ms with a time step of 0.025 ms.

.. literalinclude:: ../../python/example/single_cell_detailed_recipe.py
   :language: python
   :lines: 168-170

The results
***********

The last step is result collection. We instructed the simulation to record the spikes on the cell, and
to sample the probe.

We can print the times of the spikes:

.. literalinclude:: ../../python/example/single_cell_detailed_recipe.py
   :language: python
   :lines: 172-177

The probe results, again, warrant some more explanation:

.. literalinclude:: ../../python/example/single_cell_detailed_recipe.py
   :language: python
   :lines: 179-183

``sim.samples()`` takes a ``handle`` of the probe we wish to examine. It returns a list
of ``(data, meta)`` terms: ``data`` being the time and value series of the probed quantity; and
``meta`` being the location of the probe. The size of the returned list depends on the number of
discrete locations pointed to by the handle. We placed the probe on the "custom_terminal" locset which is
represented by 2 locations on the morphology. We therefore expect the length of ``sim.samples(handle)``
to be 2.

We plot the results using pandas and seaborn as we did in the previous example, and expect the same results:

.. literalinclude:: ../../python/example/single_cell_detailed_recipe.py
   :language: python
   :lines: 185-188

The following plot is generated. Identical to the plot of the previous example.

.. figure:: single_cell_detailed_result.svg
    :width: 400
    :align: center


The full code
*************

You can find the full code of the example at ``python/examples/single_cell_detailed_recipe.py``.
