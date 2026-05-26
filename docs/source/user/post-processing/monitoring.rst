Monitoring runs
===============

During a simulation x3d2 writes global scalar quantities to a
``monitoring.csv`` file at every output step (controlled by ``n_output``
in the input file). This file is useful for tracking the evolution of
the flow and for comparison with reference data.

Output format
-------------

The file is a comma-separated CSV with a header line:

.. code-block:: text

   # time, enstrophy, div_u_max, div_u_mean

Each subsequent row contains one record per output step.

Quantities
----------

.. list-table::
   :header-rows: 1
   :widths: 20 40 40

   * - Column
     - Description
     - Formula
   * - ``time``
     - Simulation time
     - :math:`t`
   * - ``enstrophy``
     - Spatially-averaged enstrophy
     - :math:`\mathcal{E} = \frac{1}{2N} \sum |\nabla \times \mathbf{u}|^2`
   * - ``div_u_max``
     - Maximum of :math:`|\nabla \cdot \mathbf{u}|`
     - Divergence-free check
   * - ``div_u_mean``
     - Mean of :math:`|\nabla \cdot \mathbf{u}|`
     - Divergence-free check

.. note::

   The divergence columns (``div_u_max`` and ``div_u_mean``) measure how
   well the incompressibility constraint :math:`\nabla \cdot \mathbf{u} = 0`
   is satisfied.  Since x3d2 enforces the incompressibility constraint via a pressure
   projection, these values should remain
   close to machine precision (typically :math:`\sim 10^{-14}` in double
   precision).  A growing divergence indicates a problem with the pressure
   solver or time-stepping stability.

Example: plotting with Python
-----------------------------

.. code-block:: python

   import pandas as pd
   import matplotlib.pyplot as plt

   columns = ["time", "enstrophy", "div_u_max", "div_u_mean"]
   df = pd.read_csv("monitoring.csv", comment="#", names=columns)

   fig, axes = plt.subplots(3, 1, figsize=(8, 8), sharex=True)

   axes[0].plot(df["time"], df["enstrophy"])
   axes[0].set_ylabel("Enstrophy")

   axes[1].plot(df["time"], df["div_u_max"])
   axes[1].set_ylabel("div(u) max")

   axes[2].plot(df["time"], df["div_u_mean"])
   axes[2].set_ylabel("div(u) mean")

   axes[-1].set_xlabel("Time")

   plt.tight_layout()
   plt.savefig("monitoring.png")
