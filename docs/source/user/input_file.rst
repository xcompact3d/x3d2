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
   /

``model``
  Explicit SGS model selector. Supported values are ``'none'`` and
  ``'smagorinsky'``.

  **Default:** ``'none'``

``smagorinsky_constant``
  Smagorinsky coefficient :math:`C_s`. Without wall
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

``wall_damping``
  Enables the Mason--Thomson mixing-length blend. For wall distance
  :math:`y`, roughness length :math:`z_0` and von Karman constant
  :math:`\kappa`,

  .. math::

     \ell_s = \Delta\left[C_s^{-n} +
     \left(\frac{\kappa(y+z_0)}{\Delta}\right)^{-n}\right]^{-1/n}.

  The wall distance is measured from the lower ``y`` boundary, so wall damping
  suits cases with a single wall there. The case supplies :math:`\kappa` and
  :math:`z_0`: currently only the ``abl`` case does, from its ``kappa`` and
  ``z0``, and enabling wall damping in any other case stops with an error.

  **Default:** ``.false.``

``wall_damping_n``
  Exponent :math:`n` in the Mason--Thomson blend.

  **Default:** ``3.0``

The Smagorinsky constant and ``wall_damping_n`` must be positive.

Spatial Filter
~~~~~~~~~~~~~~

An explicit low-pass filter can be applied to the velocity at the start of
every time step. It is needed for wall-modelled ABL runs: the compact schemes
cannot dissipate the :math:`2\Delta x` mode and the staggered pressure
projection cannot see it; left alone, it spoils momentum conservation in the
convection term. The filter removes it while leaving well-resolved
scales almost untouched. It is set in the ``solver_params`` namelist block.

.. code-block:: fortran

   &solver_params
     ...
     spatial_filter = .true.
     filter_alpha = 0.49
   /

``spatial_filter``
  Enables the filter. Each velocity component is filtered in
  all three directions; across a free-slip boundary the component normal to it
  is treated as odd and the other two as even.

  **Default:** ``.false.``

``filter_alpha``
  Filter parameter :math:`\alpha`, which should lie in
  :math:`(-0.5, 0.5)`. Values close to ``0.5`` act only on the smallest scales;
  smaller values damp every resolved wavenumber much harder (at ``0.40`` about
  ten times more per application than at ``0.49``). In directions that are not
  decomposed across MPI ranks the filter is solved exactly with the Thomas
  algorithm, as the system becomes only marginally diagonally dominant as
  :math:`\alpha` approaches ``0.5``.

  **Default:** ``0.49``

Atmospheric Boundary Layer Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The ``abl`` case simulates a neutral atmospheric boundary layer over a rough
wall. Select it with ``flow_case_name = 'abl'`` in ``domain_settings``, with
``y`` vertical and a no-slip floor under a free-slip lid:

.. code-block:: fortran

   &domain_settings
     flow_case_name = 'abl'
     BC_x = 'periodic', 'periodic'
     BC_y = 'dirichlet', 'neumann'
     BC_z = 'periodic', 'periodic'
     ...
   /

The case requires Smagorinsky LES with wall damping (``model = 'smagorinsky'``
and ``wall_damping = .true.`` in ``les_params``), and normally the spatial
filter. The wall damping then uses ``kappa`` and ``z0`` below. Its parameters
are set in the ``abl_nml`` namelist block.

.. code-block:: fortran

   &abl_nml
     z0 = 0.1
     u_star = 0.45
     delta = 1000.0
     kappa = 0.4
     dsampling = 3.0
     init_noise = 0.5, 0.5, 0.5
     seed = 0
     pressure_gradient = .true.
     coriolis = .false.
     coriolis_freq = 0.0
     u_geo = 0.0, 0.0, 0.0
     mass_conserve = .false.
     u_bulk = 0.0
     damping = .false.
     profile_start_iter = -1
     profile_file = 'abl_profile.csv'
   /

Driving mechanisms
^^^^^^^^^^^^^^^^^^

The flow is driven by any combination of the following switches. At least one
of ``pressure_gradient``, ``coriolis`` and ``mass_conserve`` must be on.

``pressure_gradient``
  Adds a uniform streamwise body force :math:`u_*^2/\delta`. In a steady state
  it is balanced by the wall stress.

  **Default:** ``.false.``

``coriolis``
  Adds the Coriolis force about the vertical axis,
  :math:`f\,(w,\,0,\,-u)` in :math:`(x, y, z)`. When ``pressure_gradient`` is off
  the flow is instead driven towards the geostrophic wind ``u_geo`` by the
  balancing term :math:`f\,(-U_{g,z},\,0,\,U_{g,x})`. With
  ``pressure_gradient`` on, the balancing term is left out.

  **Default:** ``.false.``

``coriolis_freq``
  Coriolis parameter :math:`f`, in 1/s.

  **Default:** ``0.0``

``u_geo``
  Geostrophic wind :math:`(U_{g,x}, U_{g,y}, U_{g,z})`. It is also the
  initial velocity when neither ``pressure_gradient`` nor ``mass_conserve`` is on,
  and the reference velocity of the damping layer for ``v`` and ``w``.

  **Default:** ``0.0, 0.0, 0.0``

``mass_conserve``
  Shifts ``u`` uniformly before every step so that the bulk velocity stays at
  its target, as the channel case does.

  **Default:** ``.false.``

``u_bulk``
  Target bulk velocity for ``mass_conserve``. If it is not positive,
  the log-law value :math:`(u_*/\kappa)\left(\ln(\delta/z_0) - \delta/L_y\right)`
  is used, where :math:`L_y` is the domain height.

  **Default:** ``0.0``

``damping``
  Adds a Rayleigh damping layer near the boundary-layer top. Its strength is
  :math:`15\,u_*/\delta`, ramping in smoothly between :math:`y = 0.95\,\delta`
  and :math:`1.05\,\delta` and applied fully above. It relaxes ``u`` towards
  :math:`(u_*/\kappa)\ln(\delta/z_0)`, and ``v`` and ``w`` towards ``u_geo``.

  **Default:** ``.false.``

Physical parameters
^^^^^^^^^^^^^^^^^^^

``u_star``
  Friction velocity :math:`u_*`, in m/s. Sets the pressure-gradient
  force, the log-law initial profile and the damping.

  **Default:** ``0.0``

``delta``
  Boundary-layer depth :math:`\delta`, in m. Must exceed ``z0``.

  **Default:** ``1.0``

``z0``
  Aerodynamic roughness length :math:`z_0`, in m. Must be positive.

  **Default:** ``0.1``

``kappa``
  von Karman constant :math:`\kappa`. Must be positive.

  **Default:** ``0.41``

Wall model
^^^^^^^^^^

The floor stress comes from the neutral log law applied at every grid column,

.. math::

   \tau_{xy} = C_d\,u\,\sqrt{u^2 + w^2}, \qquad
   \tau_{yz} = C_d\,w\,\sqrt{u^2 + w^2}, \qquad
   C_d = \left(\frac{\kappa}{\ln(h/z_0)}\right)^2,

with :math:`u` and :math:`w` sampled at height :math:`h` above the floor. It
replaces the SGS shear stress on the first grid plane above the wall, so the
wall flux passes through the same operator as the SGS stress in the interior.

``dsampling``
  Sampling height :math:`h` in grid spacings above the floor,
  rounded to the nearest grid plane. The sample must lie inside the domain and
  above ``z0``.

  **Default:** ``3.0``

Initial condition
^^^^^^^^^^^^^^^^^

The initial streamwise velocity is the log law
:math:`U(y) = (u_*/\kappa)\ln\left((y + z_0)/z_0\right)` when
``pressure_gradient`` or ``mass_conserve`` is on, otherwise the uniform
:math:`U_{g,x}`. Random noise is added, uniform in :math:`[-1, 1]` times
the ``init_noise`` amplitude: relative to :math:`U(y)` for ``u``, and in m/s
for ``v`` and ``w``. Each MPI rank draws its own random stream.

``init_noise``
  Noise amplitude for ``u``, ``v`` and ``w``. Must not be negative.

  **Default:** ``0.0, 0.0, 0.0``

``seed``
  Seed for the initial noise. With ``0`` a seed is taken from the clock,
  so every run differs. The seed used is always printed
  (``ABL noise seed: N``); setting ``seed = N`` repeats that run's initial field
  exactly with the same compiler and number of MPI ranks. Must not be negative.

  **Default:** ``0``

Profile diagnostics
^^^^^^^^^^^^^^^^^^^

The case can accumulate the horizontally and time-averaged velocity profile and
mean wall stress, sampled every ``n_output`` steps. The ``y``
direction must not be decomposed across MPI ranks.

``profile_start_iter``
  First iteration to include in the average. ``-1``
  disables the diagnostics.

  **Default:** ``-1``

``profile_file``
  CSV file rewritten at every sample with the running average.
  Its header records the sample count, the averaging interval, the imposed and
  diagnosed friction velocities (the latter from the mean wall stress) and the
  mean wall stress; its columns are ``y``, ``u_mean``, ``v_mean``, ``w_mean`` and
  the analytical log law ``u_log``.

  **Default:** ``'abl_profile.csv'``

The averages are not stored in checkpoints: after a restart the file holds the
average of the samples taken since the restart.

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
