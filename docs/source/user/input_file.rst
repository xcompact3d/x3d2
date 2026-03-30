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
     snapshot_sp = .false.
     output_stride = 2, 2, 2
     output_fields = 'pressure', 'vorticity', 'qcriterion', 'species'
     restart_from_checkpoint = .false.
     restart_file = ""
   /End

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

``snapshot_sp``: Boolean flag to save visualisation snapshots in single precision (float) instead of double precision (double). This reduces file size but may lose some precision.
**Default:** ``false`` (double precision)

``output_stride``: Three-element array specifying the spatial stride (subsampling) in ``X``, ``Y``, and ``Z`` directions for visualisation snapshots. Using values greater than ``1`` reduces file size and increases I/O performance, but decreases visualisation resolution.
  **Default:** ``[1, 1, 1]``

``output_fields``: List of additional fields to include in visualisation snapshots. Velocity components (``u``, ``v``, ``w``) are always written. Supported field names:

  - ``'pressure'`` — Pressure field, interpolated from its native cell-centred grid to the vertex grid for ParaView compatibility. Not included in checkpoint files since it is recomputed from velocity.
  - ``'vorticity'`` — Vorticity magnitude :math:`|\omega| = \sqrt{\omega_x^2 + \omega_y^2 + \omega_z^2}`, computed from the full velocity gradient tensor.
  - ``'qcriterion'`` — Q-criterion :math:`Q = -\frac{1}{2} \sum_{ij} \frac{\partial u_i}{\partial x_j} \frac{\partial u_j}{\partial x_i}`, identifying vortical structures (positive Q indicates rotation-dominated regions).
  - ``'ibm'`` — Immersed boundary method mask field (``ep1``). Values are ``1`` in the fluid domain and ``0`` in the solid domain. Requires ``ibm_on = .true.`` in the input file.
  - ``'species'`` — All transported species fields. In the input file, set ``n_species = N`` with ``N > 0`` and provide ``pr_species = ...`` in ``solver_params``, then add ``'species'`` to ``output_fields``. Snapshots then include ``phi_1`` through ``phi_N``.

  When both ``'vorticity'`` and ``'qcriterion'`` are requested, the velocity gradient tensor is computed only once.
  **Default:** (empty — only velocity is written)

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

Statistics Parameters
~~~~~~~~~~~~~~~~~~~~~

The statistics parameters control the accumulation and output of time-averaged flow statistics.
These parameters are specified in the ``stats_params`` namelist block in the input file.

.. code-block:: fortran

   &stats_params
     initstat = 0
     istatfreq = 1
     istatout = 0
     stats_prefix = "statistics"
   /

``initstat``: Timestep at which to begin accumulating statistics. Statistics are not collected before this iteration. Set to ``0`` to disable statistics entirely.
  **Default:** ``0``

``istatfreq``: Frequency (in timesteps) at which statistics are accumulated. A value of ``1`` means statistics are accumulated every timestep; a value of ``N`` means every N-th timestep.
  **Default:** ``1``

``istatout``: Frequency (in timesteps) at which accumulated statistics are written to file. Set to ``0`` to disable output.
  **Default:** ``0``

``stats_prefix``: String prefix for statistics output filenames.
  **Default:** ``"statistics"``

Output Fields
^^^^^^^^^^^^^

Statistics are written to separate ADIOS2 ``.bp`` files (e.g. ``statistics_001000.bp``), independent of the snapshot system. The output contains:

**Velocity statistics** (always present):

- ``umean``, ``vmean``, ``wmean`` — time-averaged velocity components
- ``uprime``, ``vprime``, ``wprime`` — RMS velocity fluctuations: :math:`u' = \sqrt{\max(0,\, \overline{u^2} - \bar{u}^2)}`
- ``uvmean``, ``uwmean``, ``vwmean`` — Reynolds stresses: :math:`\langle u'v' \rangle = \overline{uv} - \bar{u}\bar{v}`
- ``sample_count`` — number of samples accumulated

**Pressure statistics** (when ``'pressure'`` is in ``output_fields`` in ``checkpoint_params``):

- ``pmean`` — time-averaged pressure field

**Scalar statistics** (when ``n_species > 0`` in ``solver_params``):

- ``phimean_N`` — time-averaged scalar field for species N
- ``phiprime_N`` — RMS scalar fluctuation for species N

Accumulation uses Welford's numerically stable online algorithm. Statistics always restart from scratch — they are not saved in checkpoint files.
