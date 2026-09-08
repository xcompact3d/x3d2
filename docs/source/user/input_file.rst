Input Parameters
----------------

LES Parameters
~~~~~~~~~~~~~~

The explicit large-eddy simulation (LES) parameters are specified in the
optional ``les_params`` namelist block. Selecting ``'smagorinsky'`` computes
the eddy-viscosity field and adds the SGS stress divergence to the momentum
right-hand side at every Runge--Kutta stage or Adams--Bashforth step. This is
supported on the OpenMP and CUDA backends. If the block is omitted, LES
defaults to ``model = 'none'`` and the momentum equations are unchanged.

.. code-block:: fortran

   &les_params
     model = 'none'
     smagorinsky_constant = 0.14
     wall_damping = .false.
     wall_damping_n = 3.0
     von_karman_constant = 0.4
     roughness_length = 0.0
   /

``model``: Explicit SGS model selector. Supported values are ``'none'`` and
``'smagorinsky'``.
  **Default:** ``'none'``

``smagorinsky_constant``: Smagorinsky coefficient :math:`C_s`. Without wall
damping, the mixing length is :math:`\ell_s=C_s\Delta`, where the local filter
width is :math:`\Delta=(\Delta x\,\Delta y\,\Delta z)^{1/3}`. The eddy viscosity
is

.. math::

   \nu_t = \ell_s^2 |S|, \qquad
   |S| = \sqrt{2 S_{ij}S_{ij}}, \qquad
   S_{ij} = \frac{1}{2}\left(
     \frac{\partial u_i}{\partial x_j} +
     \frac{\partial u_j}{\partial x_i}\right).

  **Default:** ``0.14``

``wall_damping``: Enables the Mason--Thomson mixing-length blend used by
Incompact3d. For wall distance :math:`y` and roughness length :math:`z_0`,

.. math::

   \ell_s = \Delta\left[C_s^{-n} +
   \left(\frac{\kappa(y+z_0)}{\Delta}\right)^{-n}\right]^{-1/n}.

For a smooth wall (``roughness_length = 0``), this makes ``nut`` tend to zero
at the wall and approach the undamped Smagorinsky value away from it.
  **Default:** ``.false.``

``wall_damping_n``: Exponent :math:`n` in the Mason--Thomson blend.
  **Default:** ``3.0``

``von_karman_constant``: von Karman constant :math:`\kappa` used by the
wall-damping length scale.
  **Default:** ``0.4``

``roughness_length``: Aerodynamic roughness length :math:`z_0` in the same
units as the mesh coordinates. Use ``0.0`` for a smooth wall.
  **Default:** ``0.0``

The Smagorinsky, wall-damping, and von Karman constants must be positive;
``roughness_length`` must be non-negative.

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
  - ``'qcriterion'`` — Q-criterion :math:`Q = \frac{1}{2}(\nabla \cdot \mathbf{u})^2 - \frac{1}{2} \sum_{ij} \frac{\partial u_i}{\partial x_j} \frac{\partial u_j}{\partial x_i}`, identifying vortical structures (positive Q indicates rotation-dominated regions). The :math:`(\nabla \cdot \mathbf{u})^2` term is the general second-invariant form; it vanishes for divergence-free flow but is retained for numerical accuracy.
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

Wind Turbine Parameters
~~~~~~~~~~~~~~~~~~~~~~~~

The wind turbine parameters configure the ``wind_turbine`` flow case: a uniform
inflow / convective outflow domain with an optional turbine forcing model.
These parameters are specified in the ``wind_turbine_nml`` namelist block and are
only read when ``flow_case_name = 'wind_turbine'``.

.. code-block:: fortran

   &wind_turbine_nml
     init_noise = 0.0d0, 0.0d0, 0.0d0
     inlet_noise = 0.02d0, 0.0d0, 0.0d0
     bc_start_u = 8.0d0
     bc_start_v = 0.0d0
     bc_start_w = 0.0d0
     iturbine = 2
     iturboutput = 50
     adm_coords = 'discs.ad'
     rho_air = 1.2d0
     T_relax = 1.1291d0
   /End

``init_noise``: Three-element array giving the amplitude of random noise added to the whole velocity field in the initial condition, per component (``u``, ``v``, ``w``). Scaled by ``bc_start_u``.
  **Default:** ``[0, 0, 0]``

``inlet_noise``: Three-element array giving the amplitude of random noise applied to the inlet (Dirichlet) plane every substep, per component. Scaled by ``bc_start_u``. This is the recommended way to seed transition: turbulence develops physically through the wake shear layers rather than aliasing at the grid scale. A non-zero value regenerates the inlet perturbation each substep; the RNG is seeded deterministically (offset by MPI rank) so runs are reproducible for a fixed decomposition.
  **Default:** ``[0, 0, 0]``

``bc_start_u``, ``bc_start_v``, ``bc_start_w``: Uniform inflow velocity components enforced as a Dirichlet condition at the inlet plane (``x = 1``).
  **Default:** ``bc_start_u = 1``, ``bc_start_v = 0``, ``bc_start_w = 0``

``iturbine``: Selects the turbine forcing model.

  - ``0`` — no turbine forcing.
  - ``1`` — actuator line model. Not yet implemented: selecting it runs the case with no turbine forcing (equivalent to ``0``).
  - ``2`` — actuator disc model (ADM); see below.

  **Default:** ``0``

``iturboutput``: Frequency (in timesteps) at which turbine diagnostics are written. Set to ``0`` to disable.
  **Default:** ``1``

``adm_coords``: Path to the actuator-disc coordinates file (``.ad``) read when ``iturbine = 2``. See `Actuator disc file`_.
  **Default:** ``""`` (empty)

``rho_air``: Air density used in the ADM thrust and power diagnostics. The projected momentum force is divided by ``rho_air`` so it is a true acceleration; with ``rho_air = 1`` this matches the Incompact3d convention.
  **Default:** ``1``

``T_relax``: Time constant of the ADM first-order low-pass filter applied to the disc-averaged velocity. A negative value disables the filter (the instantaneous disc velocity is used).
  **Default:** ``-1`` (disabled)

Actuator Disc Model (ADM)
^^^^^^^^^^^^^^^^^^^^^^^^^^

When ``iturbine = 2``, each disc applies a thrust force computed from a thrust
coefficient and the disc-averaged inflow velocity, spread onto the grid with a
super-Gaussian smearing kernel. The thrust and power are

.. math::

   C_T' = \frac{C_T}{(1-a)^2}, \qquad
   T = \tfrac{1}{2}\, \rho_\text{air}\, C_T'\, A\, U_F^2, \qquad
   P = T\, U_F,

where :math:`A = \pi D^2/4`, :math:`a` is the induction factor, and
:math:`U_F` is the (optionally filtered) disc-averaged velocity.

Per-disc diagnostics are written to ``disc<N>.adm`` files (one per disc) with
columns for the instantaneous and filtered disc velocity, power, thrust, and
their running means. The ADM state is saved to and restored from checkpoints.

.. _Actuator disc file:

**Actuator disc file** (``.ad``): a header line followed by one row per disc:

.. code-block:: text

   CoR_x CoR_y CoR_z Yaw_deg Tilt_deg RotorDiam C_T alpha
   252.0 250.0 504.0 0.0 0.0 126.0 0.75 0.17095
   1134.0 250.0 504.0 0.0 0.0 126.0 0.75 0.17095

The columns are the disc centre (``CoR_x``, ``CoR_y``, ``CoR_z``), the yaw and
tilt angles in degrees, the rotor diameter, the thrust coefficient ``C_T``, and
the induction factor ``alpha``.

