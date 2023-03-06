Contributing
============

x3d2 is a collaborative project, open to all.  In order to enable
effective collaboration, however, we ask that your contribution(s)
comply with a set of ground rules.

In a nutshell
-------------

- For any contribution not covered by an existing issue, **open an issue
  first**.
- Respect the commit format. (Install the git hook!)
- Only commit changes formatted with `fprettify`, using the project's
  configuration.  (Install the git hook!)
- Strive to contribute **small and focused Pull Requests with detailed
  descriptions**.
- Add (a) unit test(s) covering a new feature, or (a) regression
  test(s) for a bug fix.

Issues first
------------

Issues are opened as soon as possible, particularly before any
significant implementation work is carried out or any pull request is
opened.  Starting with issues greatly facilitates collaboration:

- It encourages developers to write (and therefore think) about the
  work they plan to achieve.
- It helps maintaining a shared and up to date record of current work
  (known bugs, current developments).
- It helps coordination between developers by sign-posting who is
  working on what, as well as what is not being worked on.
- It allows discussions between developers on how to best fix the
  issue or implement the described improvement(s).

This doesn't mean that you should forbid yourself from exploratory
prototyping until an issue has been opened. In fact, prototyping is
often helpful to write the content of the issue.  On the other hand,
do open an issue before commiting to a particular design or
implementation.

Format commits accordingly
--------------------------

Commits messages are formatted according to the conventional commits
specification (v1.0.0).  Commit messages must respect the following
structure::

  <type>[optional scope]: <description>

  [optional body]

where `<type>` is one of `fix`, `feat`, `build`, `chore`, `ci`,
`docs`, `style`, `refactor`, `perf` or `test`.

Breaking changes are specified by adding an `!` before the colon `:`. Example::

  fix(allocator)!: Return 1D arrays as data blocks

  This is a breaking change because e.g. the allocator used to
  return 3D data blocks.

In addition, commit message header lines must not contain more than 68
characters.  Lines in the commit message body cannot exceed 72
characters.

.. note::

   Commit messages in x3d2 do not use footers as specified in the
   conventional commit specs.  You are welcome to use footers in your
   commit messages, but know that they hold no special place in the
   commit format policy.

A Git hook is provided to prevent you from registering non-conformant
commit messages. See setting up your development environment.

Fortran code formatting
-----------------------

