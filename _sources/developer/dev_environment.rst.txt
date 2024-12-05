.. _devenv-setup:

Setting up for development
==========================

0. To begin with, make sure you have the right tools installed by going
   through the :ref:`required tooling <tooling>`.

1. Download the x3d2 repository from GitHub::

    $ git clone git@github.com:xcompact3d/x3d2.git
    $ cd x3d2/

  The following commands assume that your shell's current directory is
  the root of the `x3d2` repository.

2. Install the `pre-commit` Git hook into your project-local
   configuration::

     $ cp githooks/pre-commit .git/hooks/
     $ chmod +x .git/hooks/pre-commit

   This Git hook will cause the automatic formatting of all fortran
   files staged in your commit, using `fprettify`.

3. Install the `commit-msg` Git hook into your project-local
   configuration::

     $ cp githooks/commit-msg .git/hooks/
     $ chmod +x .git/hooks/commit-msg

   This Git hook will automatically check your commit message against
   the commit message format, based on the conventional commits
   specification. See :ref:`the contribution guidelines
   <commit-formatting>`.
