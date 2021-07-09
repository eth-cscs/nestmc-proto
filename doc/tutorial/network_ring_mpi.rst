.. _tutorialmpi:

Distributed ring network (MPI)
==============================

In this example, the ring network created in an :ref:`earlier tutorial <tutorialnetworkring>` will be used to run the model in
a distributed context using MPI. Only the differences with that tutorial will be described.

.. Note::

   **Concepts covered in this example:**

   1. Building a basic MPI aware :py:class:`arbor.context` to run a network.
      This requires that you have built Arbor with MPI support enabled.
   2. Running the simulation and extracting the results.

The recipe
**********

Step **(11)** is changed to generate a network with five hundred cells.

.. literalinclude:: ../../python/example/network_ring_mpi.py
   :language: python
   :lines: 111-113

The hardware context
********************

Step **(12)** uses the Arbor-built-in :py:class:`MPI communicator <arbor.mpi_comm>`, which is identical to the
``MPI_COMM_WORLD`` communicator you'll know if you are familiar with MPI. The :py:class:`arbor.context` takes a
communicator for its ``mpi`` parameter. Note that you can also pass in communicators created with ``mpi4py``.
We print both the communicator and context to observe how Arbor configures their defaults.

.. literalinclude:: ../../python/example/network_ring_mpi.py
   :language: python
   :lines: 115-120

The execution
*************

Step **(16)** runs the simulation. Since we have more cells this time, which are connected in series, it will take some time for the action potential to propagate. In the :ref:`ring network <tutorialnetworkring>` we could see it takes about 5 ms for the signal to propagate through one cell, so let's set the runtime to ``5*ncells``.

.. literalinclude:: ../../python/example/network_ring_mpi.py
   :language: python
   :lines: 133-135

An important change in the execution is how the script is run. Whereas normally you run the Python script by passing
it as an argument to the ``python`` command, you need to use ``srun`` or ``mpirun`` (depending on your MPI
distribution) to execute a number of jobs in parallel. You can still execute the script using ``python``, but then
MPI will not execute on more than one node.

From the commandline, we can run the script using ``mpirun`` (``srun`` on clusters operated with SLURM) and specify the number of ranks (``NRANKS``)
or nodes. Arbor will spread the cells evenly over the ranks, so with ``NRANKS`` set to 5, we'd be spreading the 500
cells over 5 nodes, simulating 100 cells each.

.. code-block::

   mpirun -n NRANKS python mpi.py

The results
***********

Before we execute the simulation, we have to understand how Arbor distributes the computational load over the ranks.
After executing ``mpirun``, all nodes will run the same script. In the domain decomposition step, the nodes will use
the provided MPI communicator to divide the work. Once :py:func:`arbor.simulation.run` starts, each node wil work on
their allocated cell ``gid`` s.

This is relevant for the collection of results: these are not gathered for you. Remember that in step **(14)** we
store the handles to the probes; these referred to particular ``gid`` s. The ``gid`` s are now distributed, so on one
node, the script will not find the cell referred to by the handle and therefore return an empty list (no results were found).

In step **(18)** we check, for each ``gid``, if the list returned by :py:func:`arbor.simulation.samples` has a nonzero
length. The effect is that we collect the results generated on this particular node. Since we now have ``NRANKS``
instances of our script, and we can't access the results between nodes, we have to write the results to disk and
analyse them later. We query :py:attr:`arbor.context.rank` to generate a unique filename for the result.

.. literalinclude:: ../../python/example/network_ring_mpi.py
   :language: python
   :lines: 137-147

In a second script, ``network_ring_mpi_plot.py``, we load the results stored to disk into a pandas table, and plot the concatenated table as before:

.. literalinclude:: ../../python/example/network_ring_mpi_plot.py
   :language: python

To avoid an overcrowded plot, this plot was generated with just 50 cells:

.. figure:: network_ring_mpi_result.svg
    :width: 400
    :align: center

The full code
*************

You can find the full code of the example at ``python/examples/network_ring_mpi.py`` and ``python/examples/network_ring_mpi_plot.py``.
