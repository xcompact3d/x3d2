Overview of the contribution process
====================================

The graph below illustrates the general process for contributing code
changes to x3d2:

.. graphviz::

   digraph G {
       review[shape="diamond", label="Changes approved"]
       commit[label="commit changes locally"]

       PR[label="Open/update Pull Request"];

       Idea -> "Open new issue" -> commit -> PR
       PR -> review
       review -> merge [label="yes"]
       review -> commit [label=" no"]
       "Pick existing issue" -> commit

       {
           rank=sink;
           review; merge
       }
   }

Contributions are accepeted in the form of pull requests targeting the
`main` branch on the x3d2 GitHub repository.  See `(GitHub docs)
Creating a pull request
<https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request>`_.

By default your GitHub account will not be able to push changes to the
x3d2 repo, and you will have to open the pull request from a fork. See
`(GitHub docs) Creating a pull request from a fork
<https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork>`_.

Note that the whole process is **driven by issues**. If you found a
bug not currently referenced by an existing issue, or have an idea on
how to improve a part of x3d2, please open a new item on the issue
tracker before opening a pull request.
