Contributing documentation
==========================

x3d2 uses both `Sphinx <https://www.sphinx-doc.org/>`_ and `ford
<https://www.forddocs.readthedocs.io/en/latest>`_ for documentation.  Sphinx is
used to generate the user and developer guides (i.e. this website).
ford is used to generate the API documentation website by extraction
in-code documentation.

Building the documentation
--------------------------

To build the user and developer docs (Sphinx):

.. code:: console

   # From the repository root
   $ cd docs && make html

The above command generate hmtl pages in `docs/build/html`.  You can
display the user and developer docs by opening
`docs/build/html/index.html` with a web browser.


To build the API docs (ford):

.. code:: console

   # From the repository root
   $ ford ford.md -o api_docs

You can display the API docs by opening `api_docs/index.html` with a
web browser.

Writing user or developer documentation
---------------------------------------

Documentation sources are located under `docs/sources/`.  They consist
in a hierachy of reStructuredText files.  reStructuredText (rST) is a
light markup language similar to markdown.  For an introduction to
rST, see `(Sphinx)ReStructuredText primer
<https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_.

In-code documentation
---------------------

See :ref:`in-code-docs`.
