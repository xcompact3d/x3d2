.. _devenv-setup:

Setting up for development
==========================

Ensure you have the necessary tools installed by following the :ref:`required tooling <tooling>`.


1. Fork the `x3d2 repository <https://github.com/xcompact3d/x3d2/>`_ on GitHub by navigating to the x3d2 repository page and clicking the "Fork" button.

2. Clone your forked repository to your local machine:

   .. code:: bash

      $ git clone git@github.com:your-username/x3d2.git
      $ cd x3d2/

   The following commands assume that your shell's current directory is the root of the `x3d2` repository.

3. Set up the upstream repository to keep your fork in sync with the original repository:

   .. code:: bash

      $ git remote add upstream git@github.com:xcompact3d/x3d2.git
      $ git fetch upstream

4. (Optional) Install the ``pre-commit`` Git hook to automatically format Fortran files staged in your commit using ``fprettify``:

   .. code:: bash

      $ cp githooks/pre-commit .git/hooks/
      $ chmod +x .git/hooks/pre-commit