|Build| |Stars| |PyPI| |PyPIDownloads|


.. |Build| image::https://github.com/theislab/sfaira/workflows/Build%20sfaira%20Package/badge.svg
    :target: https://github.com/theislab/sfaira/workflows/Build%20sfaira%20Package/badge.svg
    :alt: Github Workflow Build sfaira Status
.. |Stars| image:: https://img.shields.io/github/stars/theislab/sfaira?logo=GitHub&color=yellow
   :target: https://github.com/theislab/sfaira/stargazers
.. |PyPI| image:: https://img.shields.io/pypi/v/sfaira?logo=PyPI
   :target: https://pypi.org/project/sfaira
.. |PyPIDownloads| image:: https://pepy.tech/badge/sfaira
   :target: https://pepy.tech/project/sfaira


sfaira - data and model repository for single-cell data
=======================================================

.. image:: https://github.com/theislab/sfaira/blob/master/resources/images/concept.png
   :width: 1000px
   :align: center

sfaira_ is a model and a data repository in a single python package (preprint_).
We provide an interactive overview of the current state of the zoos on sfaira-site_.

Its data zoo gives users access to streamlined data loaders that allow reproducible use of published and private data sets for model training and exploration.
Its model zoo gives user streamlined access to pre-trained models and to common model architectures to ease usage of neural networks in common single-cell analysis workflows:
A model zoo is a software infrastructure that improves user access to pre-trained models which are separately published, such as DCA_ or scArches_:
Instead of focussing on developing new models, we focus on making models easily accessible to users and distributable by developers.
sfaira integrates into scanpy_ workflows.

.. _scanpy: https://github.com/theislab/scanpy
.. _sfaira: https://sfaira.readthedocs.io
.. _preprint: https://www.biorxiv.org/content/10.1101/2020.12.16.419036v1
.. _DCA: https://github.com/theislab/dca
.. _scArches: https://github.com/theislab/scarches
.. _sfaira-site: https://theislab.github.io/sfaira-site/index.html
