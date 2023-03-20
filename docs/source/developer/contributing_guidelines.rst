Contributing guidelines
=======================

x3d2 is a collaborative project, open to all.  In order to enable
effective collaboration, however, we ask that your contribution(s)
comply with a set of ground rules.

In a nutshell
-------------

- For any contribution not covered by an existing issue, **open an issue
  first**.
- Respect the commit format.  (:ref:`Install the git hook! <devenv-setup>`).
- Only commit changes formatted with `fprettify`, using the project's
  configuration.  (:ref:`Install the git hook! <devenv-setup>`).
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

.. _commit-formatting:

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

Code formatting dictates things like indentation levels, number of
whitespace between operators or use upper/lowercase characters for
keywords.

Because having to learn such project specific rules isn't terribily
interesting, all Fortran code in x3d2 is formatted automatically using
`fprettify <https://github.com/pseewald/fprettify>`_. Some non-default
configuraiton options are specified in the project's `.fprettify.ini`
file.

Note that a a GitHub pull request will not be merged until all the
associated code addtions conform to the formatting.

The best way to deal with code formatting is to **automate it and
forget about it**.  To do so, make sure you install the pre-commit
hook into your local git repository, see .  This way, `fprettify` will
run automatically each time you commit new changes.

Pull requests and code review
-------------------------------

Code review is a fondamental part of x3d2's collaborative development.
Particularly, GitHub pull requests must be approved by at least two
maintainers before the associated commits are merged.  Code review allows for

- Gatekeeping, by preventing the integration of changes that are not
  up to standards, introducing bugs or potentially disruptive to
  ongoing and/or future developments.
- Knowledge transfer, by making sure that at least a couple of project
  developers have an opportunity to study and understand the code to
  be integrated.
- Training, by giving new contributors opportunities to get feedback
  and help, as well as chances to review contributions from more
  experienced developers.

Poor quality pull requests can however turn the code review process
into a bottleneck or worse, lead to poor quality reviews. In order to
facilitate the code review process, pull request must be both
**focused and small** and have a **detailed description**.

1. Pull requests should aim to merge a be coherent body of changes,
   targeting a single aspect of an issue.  It is much prefereable to
   fix an issue through multiple pull requests than one big
   contribution.
2. Pull requests should aim to merge a limited number of changes.  As
   a rule of thumb, pull requests usually become much harder to review
   beyond 10 changes.
3. Pull requests should come with a detailed description containing:

   - One or two sentences summing up the motivation for the changes.
   - A link to the addressed issue(s).
   - A summary of the changes.
   - A technical description of the changes.

   The description is important for archiving purposes, but most
   importlantl in order to guide the work of reviewers.

Note that the points above are guidelines as opposed to hard and fast
rules.  For a given pull request, the amount of work required to
review the changes will vary across reviewers.  Generally, however,
please **empathise with your fellow contributors who are going to spend
time reviewing your code**.  Aim to make their work easier, and they
will do the same for you.
