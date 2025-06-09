Contribution process
====================

The graph below illustrates the general process for contributing code
changes to x3d2:

.. mermaid::
    
   graph TD
    A[Idea] --> B["Open new issue"]
    B --> C[commit changes locally]
    C --> D[Open/update Pull Request]
    D --> E{Changes approved?}
    E -- yes --> F[merge]
    E -- no --> C
    G["Pick existing issue"] --> C

    subgraph Review
        E
        F
    end 

Contributions are accepted in the form of Pull Requests (PR) targeting the ``main`` branch on the x3d2 GitHub repository. For more information, see the GitHub documentation on `Creating a pull request <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request>`_.

By default, your GitHub account will not have push access to the x3d2 repository, so you will need to open a pull request from a fork. For details, see `Creating a pull request from a fork <https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork>`_.

Please note that the entire process is **driven by issues**. If you find a bug that is not currently referenced by an existing issue, or if you have an idea for improving x3d2, please open a new issue on the issue tracker before submitting a pull request.