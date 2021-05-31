Changelog
==========

.. role:: small
.. role:: smaller

This project adheres to `Semantic Versioning <https://semver.org/>`_.

0.2.1 :small:`2020-09-7`
~~~~~~~~~~~~~~~~~~~~~~~~

**Added**

* A commandline interface with Click, Rich and Questionary
* upgrade command, which checks whether the latest version of sfaira is installed on every sfaria startup and upgrades it if not.
* create-dataloader command which allows for the interactive creation of a sfaira dataloader script
* clean-dataloader command which cleans a with sfaira create-dataloader created dataloader script
* lint-dataloader command which runs static checks on the style and completeness of a dataloader script
* test-dataloader command which runs a unittest on a provided dataloader

**Fixed**

**Dependencies**

**Deprecated**

