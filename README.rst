|Build| |Documentation| |Stars| |PyPI| |PyPIDownloads|


.. |Build| image:: https://github.com/theislab/sfaira/workflows/Build%20sfaira%20Package/badge.svg
    :target: https://github.com/theislab/sfaira/workflows/Build%20sfaira%20Package/badge.svg
    :alt: Github Workflow Build sfaira Status

.. |Documentation| image:: https://readthedocs.org/projects/sfaira/badge/?version=latest
    :target: https://sfaira.readthedocs.io/en/latest/
    :alt: Documentation Status

.. |Stars| image:: https://img.shields.io/github/stars/theislab/sfaira?logo=GitHub&color=yellow
   :target: https://github.com/theislab/sfaira/stargazers
   :alt: Github Stars

.. |PyPI| image:: https://img.shields.io/pypi/v/sfaira?logo=PyPI
   :target: https://pypi.org/project/sfaira
   :alt: PyPI Version

.. |PyPIDownloads| image:: https://pepy.tech/badge/sfaira
   :target: https://pepy.tech/project/sfaira
   :alt: Number of downloads


sfaira - data and model repository for single-cell data
=======================================================

.. image:: https://github.com/theislab/sfaira/blob/release/resources/images/figure_rtd_intro.png
   :width: 400px
   :align: center

sfaira_ is a model and a data repository in a single python package (`full paper`_).
We provide an interactive overview of the current state of the zoos on sfaira-portal_.

Its data zoo gives users access to streamlined data loaders that allow reproducible use of published and private data sets for model training and exploration.
Its model zoo gives user streamlined access to pre-trained models and to common model architectures to ease usage of neural networks in common single-cell analysis workflows.
Instead of focussing on developing new models, we focus on making models easily accessible to users and distributable by developers.
sfaira integrates into scanpy_ workflows.

.. _scanpy: https://github.com/theislab/scanpy
.. _sfaira: https://sfaira.readthedocs.io
.. _full paper: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02452-6
.. _sfaira-portal: https://theislab.github.io/sfaira-portal/
