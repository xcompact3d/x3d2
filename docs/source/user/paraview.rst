Building ParaView with ADIOS2 for x3d2 Visualization
=====================================================

Prerequisites
-------------

- **GCC** (tested with 14.2.0)
- **CMake** (tested with 3.31.3)
- **OpenMPI** (tested with 5.0.7-GCC-14.2.0)
- **ParaView source** (tested with v6.0.1)
- **x3d2** compiled with ADIOS2 (provides ADIOS2 v2.10.2 source)
- **ffmpeg** (for video encoding)

.. note::

   x3d2 can use NVHPC + its bundled MPI for compilation. ParaView requires GCC.
   The BP file format is compiler-independent – GCC-built ParaView reads
   NVHPC-built x3d2 output without issues.

Step 1: Load Modules
--------------------

.. code-block:: bash

   module purge
   module load CMake/<version>-GCCcore-<version>
   module load OpenMPI/<version>-GCC-<version>

- Verify MPI uses GCC: ``mpicc -show`` should show ``gcc``, not ``nvc``.
- If ``OPAL_PREFIX`` is set to an NVHPC path, unset it: ``unset OPAL_PREFIX``.

Step 2: Build ADIOS2 with GCC + OpenMPI
----------------------------------------

.. code-block:: bash

   ADIOS2_SRC=<path-to-x3d2>/build-gpu-adios2/adios2-build/adios2-src
   ADIOS2_INSTALL=$HOME/adios2-gcc-install

   rm -rf $HOME/adios2-gcc-build $ADIOS2_INSTALL
   mkdir $HOME/adios2-gcc-build && cd $HOME/adios2-gcc-build

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

- Uses the ADIOS2 v2.10.2 source already downloaded by x3d2's ExternalProject.
- Installs to ``$ADIOS2_INSTALL``.

Step 3: Build ParaView
----------------------

.. code-block:: bash

   PARAVIEW_SRC=<path-to-paraview-source>
   ADIOS2_INSTALL=$HOME/adios2-gcc-install

   cd $PARAVIEW_SRC
   rm -rf build
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
     -DPARAVIEW_USE_QT=OFF \
     -DPARAVIEW_USE_PYTHON=ON \
     -DPARAVIEW_ENABLE_EXAMPLES=OFF

   make -j$(nproc)

- ``PARAVIEW_USE_QT=OFF`` builds headless (``pvpython``/``pvbatch``/``pvserver``).
- Verify cmake output shows: ``Found ADIOS2: ... (found suitable version "2.10.2")``.
- Verify ``Enabled modules: VTK(189)`` – the extra module is ``VTK::IOADIOS2``.

Key CMake notes
^^^^^^^^^^^^^^^

- ``PARAVIEW_ENABLE_ADIOS2=ON`` is **required** – without it, ``ADIOS2_DIR`` is silently ignored.
- ``ADIOS2_DIR`` must point to ``.../lib/cmake/adios2``, not just the install prefix.
- Do **not** use ``nvc++`` to compile ParaView – it fails on VTK's template-heavy C++ code.

Step 4: Verify Installation
---------------------------

.. code-block:: bash

   PVPYTHON=<path-to-paraview-source>/build/bin/pvpython

   $PVPYTHON -c 'from paraview.simple import *; print("ParaView loaded successfully")'

Test opening a BP file:

.. code-block:: bash

   cat > /tmp/test_bp.py << 'EOF'
   from paraview.simple import *
   r = OpenDataFile("<path-to-bp-file>")
   for prop in r.ListProperties():
       print(prop, "=", r.GetPropertyValue(prop))
   EOF
   $PVPYTHON /tmp/test_bp.py

Step 5: Render Movie from BP File
----------------------------------

.. code-block:: bash

   cat > /tmp/render.py << 'PYEOF'
   from paraview.simple import *
   import os

   BP_FILE = "<path-to-bp-file>"
   FIELD_NAME = "u"
   OUTPUT_DIR = "<path-to-output-frames>"

   r = ADIOS2VTXReader(FileName=BP_FILE)
   r.UpdatePipelineInformation()

   scene = GetAnimationScene()
   scene.UpdateAnimationUsingDataTimeSteps()
   tk = scene.TimeKeeper
   n_steps = len(tk.TimestepValues)
   print(f"Number of timesteps: {n_steps}")
   print(f"Time range: {tk.TimestepValues[0]} to {tk.TimestepValues[-1]}")

   view = GetActiveView() or CreateRenderView()
   view.ViewSize = [1920, 1080]

   slice1 = Slice(Input=r)
   slice1.SliceType.Normal = [0, 0, 1]
   display = Show(slice1, view)
   ColorBy(display, ('POINTS', FIELD_NAME))

   view.ResetCamera()
   Render()

   lut = GetColorTransferFunction(FIELD_NAME)
   lut.RescaleTransferFunction(-1.0, 1.0)

   os.makedirs(OUTPUT_DIR, exist_ok=True)

   SaveAnimation(os.path.join(OUTPUT_DIR, 'frame.png'), view,
       ImageResolution=[1920, 1080],
       FrameWindow=[0, n_steps - 1])

   print(f"Done. {n_steps} frames saved to {OUTPUT_DIR}")
   PYEOF
   $PVPYTHON /tmp/render.py

Step 6: Encode to MP4
---------------------

.. code-block:: bash

   ffmpeg -framerate 24 -i <path-to-output-frames>/frame.%04d.png \
     -c:v libx264 -profile:v high -crf 20 -pix_fmt yuv420p \
     <path-to-output>/movie.mp4

- ``framerate 24`` → ~20s video for 477 frames. Use ``60`` for faster playback.
- ``crf 20`` → good quality. Lower = better quality, larger file.
- Download the mp4: ``scp user@host:<path-to-output>/movie.mp4 .``

Troubleshooting
---------------

- **ADIOS2_DIR unused warning:** You forgot ``-DPARAVIEW_ENABLE_ADIOS2=ON``.
- **MPI init errors with /proj/nv/... paths:** ``OPAL_PREFIX`` is set to NVHPC path. Run ``unset OPAL_PREFIX``.
- **nvc++ compile errors in VTK:** Use GCC instead. Rebuild ADIOS2 with GCC + OpenMPI.
- **pvpython hangs:** Check MPI environment. Try ``export OMPI_MCA_btl=self,tcp``.
- **No timesteps / wrong FrameWindow:** Use ``len(tk.TimestepValues)`` to get the actual count.
