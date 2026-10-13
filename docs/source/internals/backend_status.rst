Backend implementation status
=============================

Every operation the solver needs is a type-bound procedure of
``base_backend_t`` (``src/backend/backend.f90``). None of them is
deferred: the base class provides a default that stops the run with
``<operation>: not implemented by the active backend``, and each backend
overrides only the operations it actually implements.

The table below lists, for each operation, whether the OpenMP CPU
(``omp``), CUDA (``cuda``) and OpenMP target offload (``omptgt``)
backends override it.

.. list-table::
   :header-rows: 1
   :widths: 30 10 10 10

   * - Operation
     - OMP
     - CUDA
     - OMP_TGT
   * - ``transeq_x``
     - ✅
     - ✅
     - ✅
   * - ``transeq_y``
     - ✅
     - ✅
     - ✅
   * - ``transeq_z``
     - ✅
     - ✅
     - ✅
   * - ``transeq_species``
     - ✅
     - ✅
     - ✅
   * - ``tds_solve``
     - ✅
     - ✅
     - ✅
   * - ``thom_solve``
     - ✅
     - ✅
     - ✅
   * - ``reorder``
     - ✅
     - ✅
     - ✅
   * - ``sum_yintox``
     - ✅
     - ✅
     - ✅
   * - ``sum_zintox``
     - ✅
     - ✅
     - ✅
   * - ``veccopy``
     - ✅
     - ✅
     - ✅
   * - ``vecadd``
     - ✅
     - ✅
     - ✅
   * - ``vecmult``
     - ✅
     - ✅
     - ✅
   * - ``scalar_product``
     - ✅
     - ✅
     - ✅
   * - ``vector_norm_squared``
     - ✅
     - ✅
     - ✅
   * - ``field_max_mean``
     - ✅
     - ✅
     - ✅
   * - ``slice_max_sum``
     - ✅
     - ✅
     - ✅
   * - ``field_plane_sums``
     - ✅
     - ✅
     - ✅
   * - ``field_scale``
     - ✅
     - ✅
     - ✅
   * - ``field_shift``
     - ✅
     - ✅
     - ✅
   * - ``field_volume_integral``
     - ✅
     - ✅
     - ✅
   * - ``field_set_face``
     - ✅
     - ✅
     - ✅
   * - ``field_set_y_plane``
     - ✅
     - ✅
     - ✅
   * - ``field_set_abl_wall_stress``
     - ✅
     - ✅
     - ✅
   * - ``field_set_face_from_field``
     - ✅
     - ✅
     - ✅
   * - ``field_add_face_from_field``
     - ✅
     - ✅
     - ✅
   * - ``compute_vorticity``
     - ✅
     - ✅
     - ✅
   * - ``compute_qcriterion``
     - ✅
     - ✅
     - ✅
   * - ``compute_smagorinsky_nut``
     - ✅
     - ✅
     - ✅
   * - ``compute_sgs_stress``
     - ✅
     - ✅
     - ✅
   * - ``copy_data_to_f``
     - ✅
     - ✅
     - ✅
   * - ``copy_f_to_data``
     - ✅
     - ✅
     - ✅
   * - ``alloc_tdsops``
     - ✅
     - ✅
     - ✅
   * - ``init_poisson_fft``
     - ✅
     - ✅
     - ❌ [#fft]_
   * - ``sync``
     - ✅
     - ✅
     - ✅ [#sync]_
   * - ``get_device_bw_info``
     - ✅
     - ✅
     - ✅ [#bw]_
   * - ``supports_device_field_export``
     - ➖ [#export]_
     - ✅
     - ➖ [#export]_
   * - ``export_field_to_device``
     - ➖ [#export]_
     - ✅
     - ➖ [#export]_

Legend: ✅ implemented, ❌ not implemented (falls back to the
``not_implemented`` default), ➖ base-class default is the correct
behaviour.

``base_init``, ``get_field_data`` and ``set_field_data`` are concrete in
``base_backend_t`` and shared by all backends; they work wherever
``reorder``, ``copy_data_to_f`` and ``copy_f_to_data`` are implemented.

.. [#sync] Runs an empty target region that maps a result back, which the
   host must wait for: some runtimes (CCE) return from a target region
   before its kernel completes, and OpenMP has no direct call to wait for
   the device.
.. [#bw] OpenMP has no query for the device memory clock or bus width, so
   ``omptgt`` reads them from the vendor runtime (CUDA driver on NVIDIA, HIP
   on AMD), opened at run time with ``dlopen``; it reports
   ``available = .false.`` when that runtime is not found. See
   ``src/backend/omptgt/device_info.f90``.
.. [#fft] Not yet implemented: the CUDA backend builds its Poisson solver
   on cuFFT/cuFFTMp, and ``omptgt`` has no FFT library behind it yet. Until
   it does, ``omptgt`` cannot run a flow case.
.. [#export] Not overridden. The default ``supports_device_field_export``
   returns ``.false.``, so ``export_field_to_device`` is never called.
