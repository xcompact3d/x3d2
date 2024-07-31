.. _tooling:

Required tools
==============

To contribute to x3d2 you need a few things on your toolbelt

- A Fortran 2003 compiler (e.g. gfortran)
- The `NVIDA HPC Fortran compiler <https://developer.nvidia.com/hpc-compilers>`_ [optional]
- The CMake build system
- The fprettify auto-formatter
- The FORD and Sphinx documentation generators.

The above tools can all be installed via you distribution's package
manager or the Python package pamanger `pip`.

.. note::

   We strongly recommend to use the `pipx` wrapper around `pip` to
   ensure installed Python packages are installed into isolated
   virtual environments.

Fortran compiler
----------------

You can get started with the GNU Fortran compiler

.. code:: console

   $ sudo apt install gfortran

CMake
-----

CMake is used to configure the build.  This means, among other things,
finding the location of the required compiler and libraries, as well
as resolving the links between different components of the software.

To configure and build x3d2 you will need CMake version 3.18 and
above.  At the time of writing the currently distributed version of
CMake in Debain stable is 3.18, so you should be able to install CMake
from your distribution's package manager, for instance

.. code:: console

   sudo apt install cmake

Recent CMake versions are packaged for a variety of package managers,
including `pipx`:

.. code:: console

   $ pipx install cmake

fprettify
---------

Auto-formatter for fortran 90 and above.

.. code:: console

   $ pipx install fprettify

FORD
----

Documentation generator for fortran 90 and above.

.. code:: console

   $ pipx install ford

Sphinx
------

.. code:: console

   $ pipx install sphinx


