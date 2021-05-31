Installation
============

sfaira is pip installable.

PyPI
~~~~
To install a sfaira release directly from PyPi, run::

    pip install sfaira


Install a development version
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To install a specific branch `target_branch` of sfaira from a clone, run::

    cd target_directory
    git clone https://github.com/theislab/sfaira.git
    cd sfaira
    git checkout target_branch
    git pull
    pip install -e .

In most cases, you would install one of the following:
You may choose the branch `release` if you want to use a relatively stable version
which is similar to the current release but may have additional features already.
You may choose the branch `dev` if you want newer features than available from `release`.
You may choose a specific feature branch if you want to use or improve that feature before it
is reviewed and merged into `dev`.
