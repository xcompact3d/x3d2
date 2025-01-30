Contributing guidelines
=======================

We welcome your contributions to x3d2! To ensure effective collaboration, please follow the guidelines outlined below.


General principles
~~~~~~~~~~~~~~~~~~

* Create an issue before opening a Pull Request (PR)
* Keep PRs small and focused with detailed descriptions
* Ensure Fortran code is formatted with ``fprettify``
* Write tests for new features and bug fixes
* Request a maintainer to review your PR before merging
* Rebase your branch on top of the latest ``main`` before merging


Issues first
~~~~~~~~~~~~

* Open an issue before starting work on a feature or bug fix, unless one already exists.
* Open issues as early as possible to avoid wasted effort and ensure alignment.
* Provide detailed descriptions and relevant context in issues.
* Assign issues to the person responsible for the work to avoid duplication.
* Keep issues updated with the latest information and status.
* Close issues once the work is complete to maintain project organization.

.. note::

  Prototyping is encouraged and can help write the issue content. However, open an issue before committing to a specific design or implementation.

.. _commit-formatting:

Commit message guidelines
~~~~~~~~~~~~~~~~~~~~~~~~~

To keep commit history readable and maintainable, we recommend the following best practices:

* Write clear and concise commit messages.
* Use short messages for minor changes. Messages like "Fix typo" or "Update documentation" are acceptable for small updates.
* Make atomic commits: Each commit should represent a single logical change. Avoid bundling unrelated changes in a single commit.
* Reference relevant issue numbers in commit messages to provide context (e.g., "Fixes #123").
* Provide a detailed description in the commit message body if the change is not self-explanatory.


Code formatting
~~~~~~~~~~~~~~~

* Use `fprettify <https://github.com/pseewald/fprettify>`_ for automatic Fortran formatting. Some non-default configuraiton options are specified in the project's ``.fprettify.ini`` file. 
* We recommend installing the pre-commit hook to automatically format your code before each commit. Alternatively you can run ``fprettify`` manually before submitting a pull request.
* PRs will not be merged unless all code conforms to formatting.

Pull requests and code review
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


* PRs must be small and focused (one PR per feature/bug fix).
* Write meaningful PR desriptions to help reviewers. This should include:

  * Summary of changes
  * Technical explanation of key modifications
  * Link to the relevant issue
  * Test results (if applicable)

* Self-review your code and ensure tests pass before requesting a review.
* PRs must be approved by at least one maintainer before merging.
* For high-impact changes a second review is recommended but not required.
* For small, non-critical PRs (e.g. documentation, minor formatting), self-merging is allowed.

Merging best practices
~~~~~~~~~~~~~~~~~~~~~~

Before merging your PR to ``main``,  follow these guidelines.

* Rebase instead of merging: Rebase your feature branch on top of the latest ``main`` to avoid unnecessary merge commits and keep the commit history linear.

.. code-block:: bash

    # Fetch the latest changes from the original repository
    git fetch upstream

    # Checkout your feature branch
    git checkout my-feature-branch

    # Rebase your branch on top of the latest main
    git rebase upstream/main

* Resolve conflicts during rebase: If there are conflicts during rebase, resolve them manually. Use the following commands to continue the rebase after resolving conflicts:

.. code-block:: bash

    # Add the resolved files
    git add .

    # Continue the rebase
    git rebase --continue

* Squash commits before merging: If your feature branch has multiple commits, squash them into a single commit before merging. This keeps the commit history clean and concise.