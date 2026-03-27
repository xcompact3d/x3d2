Building ParaView with ADIOS2
=============================

x3d2 writes output in the `ADIOS2 BP format <https://adios2.readthedocs.io/>`_.
To visualise these files in `ParaView <https://www.paraview.org/>`_, ParaView must
be built from source with ADIOS2 support enabled. This guide walks through
building ADIOS2 and ParaView with the flags required to read x3d2 output.

.. note::

   If you already have a ParaView installation with ADIOS2 support (check
   *Help → About* for the ``VTK::IOADIOS2`` module), you can skip this guide
   entirely.

Prerequisites
-------------

You will need the following tools installed on your system:

- A C/C++ compiler — `GCC <https://gcc.gnu.org/install/>`_ is recommended
  (tested with GCC 14.2). Do not use ``nvc++``; it fails on VTK's
  template-heavy C++ code.
- `CMake <https://cmake.org/download/>`_ 3.12 or newer (tested with 3.31.3).
- MPI — any MPI implementation that uses your C/C++ compiler
  (e.g. `OpenMPI <https://docs.open-mpi.org/en/v5.0.x/installing-open-mpi/quickstart.html>`_,
  MPICH). Verify with ``mpicc -show``; it should report
  ``gcc``, not ``nvc``.
- Python 3 (3.3+) — required if you want to script ParaView with
  ``pvpython``.
- Qt 5.15+ — required only if you want the ParaView GUI. Omit if you only
  need CLI scripting (``pvpython`` / ``pvbatch``).

.. note::

   x3d2 can be compiled with NVHPC and its bundled MPI. ParaView requires
   GCC. This is not a problem — the BP file format is compiler-independent,
   so a GCC-built ParaView reads NVHPC-built x3d2 output without issues.

.. tip::

   On HPC systems that use environment modules, load GCC-based toolchains
   before building (e.g. ``module load CMake OpenMPI``). Make sure
   ``OPAL_PREFIX`` is not set to an NVHPC path — unset it if necessary
   (``unset OPAL_PREFIX``).

Step 1: Build ADIOS2
---------------------

Full installation instructions are available in the
`ADIOS2 installation guide <https://adios2.readthedocs.io/en/latest/setting_up/setting_up.html>`_.
Below is a summary of the steps with the specific options needed for use with
ParaView and x3d2.

Clone the ADIOS2 source:

.. code-block:: bash

   git clone https://github.com/ornladios/ADIOS2.git ADIOS2

Configure and build:

.. code-block:: bash

   ADIOS2_SRC=$PWD/ADIOS2
   ADIOS2_INSTALL=$HOME/adios2-install

   mkdir adios2-build && cd adios2-build

   cmake $ADIOS2_SRC \
     -DCMAKE_C_COMPILER=gcc \
     -DCMAKE_CXX_COMPILER=g++ \
     -DCMAKE_Fortran_COMPILER=mpifort \
     -DMPI_C_COMPILER=mpicc \
     -DMPI_CXX_COMPILER=mpicxx \
     -DADIOS2_USE_Fortran=ON \
     -DADIOS2_USE_CUDA=OFF \
     -DADIOS2_USE_MPI=ON \
     -DADIOS2_BUILD_EXAMPLES=OFF \
     -DADIOS2_BUILD_TESTING=OFF \
     -DCMAKE_INSTALL_PREFIX=$ADIOS2_INSTALL

   make -j$(nproc)
   make install

The key options here are:

- ``ADIOS2_USE_Fortran=ON`` — needed because x3d2 is a Fortran code.
- ``ADIOS2_USE_CUDA=OFF`` — CUDA is not needed for the ParaView reader.
- ``ADIOS2_USE_MPI=ON`` — enables parallel I/O.

Step 2: Build ParaView
-----------------------

Full build instructions for all platforms (Linux, macOS, Windows) are available
in the `ParaView Build Guide <https://www.paraview.org/paraview-docs/latest/cxx/md__builds_gitlab-kitware-sciviz-ci_Documentation_dev_build.html>`_.
Below are the specific CMake flags needed to enable ADIOS2 support for x3d2.

Clone ParaView (replace ``v6.0.1`` with your desired version):

.. code-block:: bash

   git clone https://gitlab.kitware.com/paraview/paraview.git
   cd paraview
   git checkout v6.0.1
   git submodule update --init --recursive

Configure and build:

.. code-block:: bash

   ADIOS2_INSTALL=$HOME/adios2-install

   mkdir build && cd build

   cmake .. \
     -DCMAKE_BUILD_TYPE=Release \
     -DCMAKE_C_COMPILER=gcc \
     -DCMAKE_CXX_COMPILER=g++ \
     -DPARAVIEW_USE_MPI=ON \
     -DMPI_C_COMPILER=mpicc \
     -DMPI_CXX_COMPILER=mpicxx \
     -DPARAVIEW_ENABLE_ADIOS2=ON \
     -DADIOS2_DIR=$ADIOS2_INSTALL/lib/cmake/adios2 \
     -DPARAVIEW_USE_PYTHON=ON

   make -j$(nproc)

.. important::

   - ``PARAVIEW_ENABLE_ADIOS2=ON`` is required. Without it,
     ``ADIOS2_DIR`` is silently ignored.
   - ``ADIOS2_DIR`` must point to ``<install-prefix>/lib/cmake/adios2``, not
     just the install prefix.

With or without the GUI:

- To build the full ParaView GUI, install Qt 5.15+ and ensure it is
  discoverable by CMake (the default ``PARAVIEW_USE_QT=ON`` will find it).
- To build CLI (scripting only via ``pvpython`` / ``pvbatch`` /
  ``pvserver``), add ``-DPARAVIEW_USE_QT=OFF``. This is useful on HPC
  clusters without display servers.

Verifying the build:

After CMake configuration completes, check the output for:

- ``Found ADIOS2: ... (found suitable version ...)`` — confirms ADIOS2 was
  detected.
- ``VTK::IOADIOS2`` listed among the enabled VTK modules.

Step 3: Verify Installation
----------------------------

Test that ParaView loads correctly:

.. code-block:: bash

   ./build/bin/pvpython -c 'from paraview.simple import *; print("ParaView OK")'

Test opening a BP file produced by x3d2:

.. code-block:: bash

   cat > /tmp/test_bp.py << 'EOF'
   from paraview.simple import *
   r = OpenDataFile("/path/to/your/file.bp")
   for prop in r.ListProperties():
       print(prop, "=", r.GetPropertyValue(prop))
   EOF
   ./build/bin/pvpython /tmp/test_bp.py

Replace ``/path/to/your/file.bp`` with the path to an actual BP file from an
x3d2 simulation. If the script prints the file's properties without errors,
the installation is working.

Troubleshooting
---------------

- ``ADIOS2_DIR`` unused warning — You forgot ``-DPARAVIEW_ENABLE_ADIOS2=ON``.
- MPI initialisation errors mentioning ``/proj/nv/...`` — ``OPAL_PREFIX`` is
  set to an NVHPC path. Run ``unset OPAL_PREFIX``.
- ``nvc++`` compile errors in VTK — Use GCC. Rebuild both ADIOS2 and
  ParaView with GCC.
- ``pvpython`` hangs — Check your MPI environment. Try
  ``export OMPI_MCA_btl=self,tcp``.