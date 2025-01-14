Getting Started
===============

Dependencies
------------

NVIDIA HPC SDK
~~~~~~~~~~~~~~

To install x3d2 for NVIDIA GPUs - download and install the NVIDIA HPC SDK for your target platform from the `NVIDIA website <https://developer.nvidia.com/hpc-sdk-downloads>`_

CMake
~~~~~

To build x3d2, you will need CMake version 3.20 and above. You can download the latest version from the `CMake website <https://cmake.org/download/>`_

Installing on Linux
-------------------

To install x3d2 from source, follow these steps:

1. Clone the repository and change to x3d2 directory

.. code-block:: bash

   git clone https://github.com/xcompact3d/x3d2.git
   cd x3d2

2. Set the path variables as follows (consider adding them to your ``~/.bashrc`` file if you want them to persist)

.. code-block:: bash

   export PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.11/comm_libs/mpi/bin:/opt/nvidia/hpc_sdk/Linux_x86_64/24.11/compilers/bin:$PATH
   export LD_LIBRARY_PATH=/opt/nvidia/hpc_sdk/Linux_x86_64/24.11/comm_libs/mpi/lib:/opt/nvidia/hpc_sdk/Linux_x86_64/24.11/compilers/lib:$LD_LIBRARY_PATH

3. Ensure that ``mpirun`` and ``mpif90`` picks up the version shipped with NVIDIA HPC SDK. For instance it should give:

.. code-block:: bash

   $ which mpirun
   /opt/nvidia/hpc_sdk/Linux_x86_64/24.11/comm_libs/mpi/bin/mpirun

   $ which mpif90
   /opt/nvidia/hpc_sdk/Linux_x86_64/24.11/comm_libs/mpi/bin/mpif90

4. Set Fortran Compiler flag to use ``mpif90``

.. code-block:: bash

   export FC=mpif90

5. Create the build system using CMake

.. code-block:: bash

   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

where ``build`` specifies the directory to which we write the build configuration files and ``DCMAKE_BUILD_TYPE=Release`` specifies the build type (in this case Release build; if you want to install the debug please use ``Debug`` instead).

6. Change into the build directory and create the build

.. code-block:: bash

   cd build
   make

This should create a binary file called ``xcompact`` within ``build/src/`` directory.