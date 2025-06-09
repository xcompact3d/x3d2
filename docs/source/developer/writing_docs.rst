Contributing to documentation
=============================

x3d2 uses both `Sphinx <https://www.sphinx-doc.org/>`_ and `FORD <https://www.forddocs.readthedocs.io/en/latest>`_ for documentation. 

- **Sphinx** is used to generate the user and developer guides (i.e., this website).
- **FORD** is used to generate the API documentation by extracting in-code comments.

To contribute to the documentation, ensure you have both tools installed. See :ref:`tooling` for installation instructions.

Building documentation
----------------------

To build the user and developer guides (Sphinx):

.. code:: bash

   # From the repository root
   $ cd docs && make html

The above command generates HTML pages in ``docs/build/html``. You can view the user and developer guides by opening ``docs/build/html/index.html`` in a web browser.

To build the API documentation (FORD):

.. code:: bash

   # From the repository root
   $ ford ford.md -o api_docs

You can view the API documentation by opening ``api_docs/index.html`` in a web browser.

Writing Documentation
---------------------

Documentation sources are located under ``docs/source/``. They consist of a hierarchy of reStructuredText (rST) files. rST is a lightweight markup language similar to Markdown. For an introduction to rST, see the `(Sphinx) reStructuredText Primer <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_.

In-code Documentation
---------------------

For information on writing in-code documentation, see :ref:`in-code-docs`.