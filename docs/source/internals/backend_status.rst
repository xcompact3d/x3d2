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
     - ❌
   * - ``transeq_y``
     - ✅
     - ✅
     - ❌
   * - ``transeq_z``
     - ✅
     - ✅
     - ❌
   * - ``transeq_species``
     - ✅
     - ✅
     - ❌
   * - ``tds_solve``
     - ✅
     - ✅
     - ❌
   * - ``thom_solve``
     - ✅
     - ✅
     - ❌
   * - ``reorder``
     - ✅
     - ✅
     - ✅
   * - ``sum_yintox``
     - ✅
     - ✅
     - ❌
   * - ``sum_zintox``
     - ✅
     - ✅
     - ❌
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
     - ❌
   * - ``scalar_product``
     - ✅
     - ✅
     - ❌
   * - ``vector_norm_squared``
     - ✅
     - ✅
     - ✅
   * - ``field_max_mean``
     - ✅
     - ✅
     - ❌
   * - ``slice_max_sum``
     - ✅
     - ✅
     - ❌
   * - ``field_scale``
     - ✅
     - ✅
     - ❌
   * - ``field_shift``
     - ✅
     - ✅
     - ❌
   * - ``field_volume_integral``
     - ✅
     - ✅
     - ❌
   * - ``field_set_face``
     - ✅
     - ✅
     - ❌
   * - ``field_set_face_from_field``
     - ✅
     - ✅
     - ❌
   * - ``compute_vorticity``
     - ✅
     - ✅
     - ❌
   * - ``compute_qcriterion``
     - ✅
     - ✅
     - ❌
   * - ``compute_smagorinsky_nut``
     - ✅
     - ✅
     - ❌
   * - ``compute_sgs_stress``
     - ✅
     - ✅
     - ❌
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
     - ❌
   * - ``init_poisson_fft``
     - ✅
     - ✅
     - ❌
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

.. [#sync] No-op: every target region in ``omptgt`` is synchronous.
.. [#bw] Always reports ``available = .false.``: OpenMP has no portable
   query for the device memory clock or bus width.
.. [#export] Not overridden. The default ``supports_device_field_export``
   returns ``.false.``, so ``export_field_to_device`` is never called.
