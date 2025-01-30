.. _tooling:

Required tools
==============

To contribute to x3d2, ensure you have the following tools:

- A Fortran 2003 compiler (e.g. gfortran)
- The `NVIDA HPC SDK  <https://developer.nvidia.com/hpc-sdk>`_ (optional)
- `CMake <https://cmake.org/>`_ build system
- `fprettify <https://github.com/pseewald/fprettify>`_ auto-formatter
- Open MPI (if not using NVIDA HPC SDK)
- The `FORD <https://github.com/Fortran-FOSS-Programmers/ford>`_ and `Sphinx <https://www.sphinx-doc.org/en/master/>`_ documentation generators.
- `Git <https://git-scm.com/>`_ for version control

The above tools can be installed via your platform's package manager. Note that the NVIDIA HPC SDK must be downloaded and installed from the NVIDIA website.

.. note::

   We recommend using the ``pipx`` wrapper around ``pip`` to ensure Python packages are installed into isolated virtual environments.

Fortran compiler
----------------

You can get started with the GNU Fortran compiler. On Ubuntu, you can

.. code:: console

   $ sudo apt install gfortran

CMake
-----

CMake is used to configure the build process, including locating the required compilers and libraries, and resolving dependencies between different components of the software.

To configure and build x3d2, you will need CMake version 3.18 or above. You should be able to install CMake from your distribution's package manager. For example, on Debian-based systems, you can install CMake with:

.. code:: console

   $ sudo apt install cmake

Alternatively, you can install CMake using ``pipx``:

.. code:: console

   $ pipx install cmake

Auto-formatter for Fortran 90 and above
---------------------------------------

To ensure consistent formatting of Fortran code, we use `fprettify`. You can install it using ``pipx``:

.. code:: console

   $ pipx install fprettify

API documentation generator
---------------------------

To generate documentation for Fortran code, we use FORD. You can install it using ``pipx``:


.. code:: console

   $ pipx install ford


General documentation generator
-------------------------------

To generate general documentation, we use `Sphinx`. You can install it using ``pipx``:

.. code:: console

   $ pipx install sphinx