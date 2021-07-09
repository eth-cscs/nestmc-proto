.. _contribdoc:

Documentation
=============

The source for `the documentation <https://docs.arbor-sim.org>`__ is
found in the ``/doc`` subdirectory. You can add your contribution to the documentation
in the same way you would contribute code, please see the :ref:`contribpr` section.

.. _contribdoc-tut:

Tutorials
---------

Documentation includes tutorials, which are :ref:`examples <contribexample>` with an associated
page under ``/doc/tutorial`` explaining the goal of the example and guiding the reader through the code.
We try to use the `literalinclude directive <https://www.sphinx-doc.org/en/master/usage/restructuredtext/directives.html#directive-literalinclude>`
wherever possible to ensure the tutorial and example don't diverge.
Remember to update the line numbers whenever you update the examples in the tutorial.

Update policy
-------------

How to we decide if documentation is good? By observing how effective it is used
in practice. If a question on Arbor (regardless of medium) is satisfactorily
resolved (on both sides!) by with a pointer to the (relevant section in the) docs,
the documentation was good. If, on the other hand, explanation was needed, the
documentation was bad.

If you found the documentation to be insufficiently clear or elaborate, you must
consider this a bug and find or file an `issue <https://github.com/arbor-sim/arbor/issues>`__ and if you are able, make a :ref:`pull request <contribpr>`.

.. _contribdoc-namingconventions:

Labels and filename conventions
-------------------------------

Although it is not absolutely essential to do so, we try to keep to the following conventions
for naming files and labels, with the goal of making it easy to construct one from the other
such that you don't have to remember or look anything up. We try to cross-link where we can,
which means we label where we can, which translates to a large number of labels.

Wherever possible, observe:

* file names: avoid underscores as much as possible. E.g. `cpp/cable_cell.rst` -> `cpp/cablecell.rst`.
* page-level labels in each file: the path concatenated without spaces. E.g. `cpp/cablecell.rst` -> `cppcablecell`.
* heading labels in a file: `pagelevel-sectionname`. Feel free to slugify long headings.
E.g. the morphology section in the C++ cable cell docs could be `cppcablecell-morph`.

.. _contribdoc-lang:

Language
--------

Although the primary language of Arbor is C++, we use English for most of the non-code.
If in doubt, we recommend following the European Union's
`English style guide <https://ec.europa.eu/info/sites/info/files/styleguide_english_dgt_en.pdf>`_.

In general we try to have a relaxed and concise approach to the language.
