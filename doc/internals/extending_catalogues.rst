.. _extending-catalogues:

Adding Built-in Catalogues to Arbor
===================================

There are two ways new mechanisms catalogues can be added to Arbor, statically
or dynamically. If you have a set of mechanisms to use with Arbor, you are in
all likelihood interested in the former.

.. warning::

   If you are coming from NEURON and looking for the equivalent of
   ``nrnivmodl``, please read on :ref:`here <_mechanisms_dynamic>`.

   Following this path is for developers rather than end-users.

This requires a copy of the Arbor source tree and the compiler toolchain used to
build Arbor in addition to the installed library. Following these steps will
produce a catalogue of the same level of integration as the built-in catalogues
(*default*, *bbp*, and *allen*). The required steps are as follows

1. Go to the Arbor source tree.
2. Create a new directory under *mechanisms*.
3. Add your .mod files.
4. Edit *mechanisms/CMakeLists.txt* to add a definition like this (example from
   *default* catalogue)

   .. code-block :: cmake

     make_catalogue(
       NAME default                                                   # Name of your catalogue
       SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/default"                  # Directory name (added above)
       OUTPUT "CAT_DEFAULT_SOURCES"                                   # Variable name to store C++ files into (see below)
       MECHS exp2syn expsyn expsyn_stdp hh kamt kdrmt nax nernst pas  # Space separated list of mechanisms w/o .mod suffix.
       PREFIX "${PROJECT_SOURCE_DIR}/mechanisms"                      # where does 'generate_catalogue' live, do not change
       STANDALONE FALSE                                               # build as shared object, must be OFF
       VERBOSE OFF)                                                   # Print debug info at configuration time

5. Add your ``output-name`` to the ``arbor_mechanism_sources`` list.

   .. code-block :: cmake

     set(arbor_mechanism_sources
       ${CAT_BBP_SOURCES}
       ${CAT_ALLEN_SOURCES}
       ${CAT_DEFAULT_SOURCES}                                          # from above
       PARENT_SCOPE)

6. Add a ``global_NAME_catalogue`` function in ``mechcat.hpp`` and ``mechcat.cpp``
7. Bind this function in ``python/mechanisms.cpp``.

All steps can be more or less copied from the surrounding code.
