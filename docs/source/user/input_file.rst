Input Parameters
----------------

Checkpoint Parameters
~~~~~~~~~~~~~~~~~~~~~

The checkpoint parameters control how the simulation saves its state for both restart purposes and visualisation. 
These parameters are specified in the ``checkpoint_params`` namelist block in the input file.

.. code-block:: fortran

   &checkpoint_params
     checkpoint_freq = 1000
     snapshot_freq = 500
     keep_checkpoint = .false.
     checkpoint_prefix = "checkpoint"
     snapshot_prefix = "snapshot"
     output_stride = 2, 2, 2
     restart_from_checkpoint = .false.
     restart_file = ""
   /

``checkpoint_freq``: Frequency (in timesteps) at which to save checkpoint files for simulation restart. Set to ``0`` to disable checkpointing.
  **Default:** ``0``

``snapshot_freq``: Frequency (in timesteps) at which to save visualisation snapshot files. Set to ``0`` to disable snapshots.
  **Default:** ``0``

``keep_checkpoint``: Controls whether to keep all checkpoint files (``true``) or only the most recent one (``false``).
  **Default:** ``false``

``checkpoint_prefix``: String prefix for checkpoint filenames. Each checkpoint will be named as ``<checkpoint_prefix>_XXXXXX.bp`` where ``XXXXXX`` is the timestep number.
  **Default:** ``"checkpoint"``

``snapshot_prefix``: String prefix for visualisation snapshot filenames. Each snapshot will be named as ``<snapshot_prefix>_XXXXXX.bp``.
  **Default:** ``"snapshot"``

``output_stride``: Three-element array specifying the spatial stride (subsampling) in ``X``, ``Y``, and ``Z`` directions for visualisation snapshots. Using values greater than ``1`` reduces file size and increases I/O performance, but decreases visualisation resolution.
  **Default:** ``[1, 1, 1]``

``restart_from_checkpoint``: Boolean flag to restart the simulation from a checkpoint file.
  **Default:** ``false``

``restart_file``: Path to the checkpoint file for restarting the simulation. Required when ``restart_from_checkpoint`` is ``true``.
  **Default:** ``""`` (empty string)

Technical Details
^^^^^^^^^^^^^^^^^

The checkpoint system uses ADIOS2 BP format for I/O which provides:

- Efficient parallel I/O even on large HPC systems
- Compression options to reduce storage requirements
- Compatibility with visualisation tools (ParaView can directly read BP files)
- Restart files contain full-resolution field data
- Visualisation files can use strided (lower resolution) output for performance

To view snapshot files in ParaView, open the generated ``.bp`` files using the ADIOS2 reader plugin.
