# -*- coding: utf-8 -*-
"""A Data and Model Zoo for Single-Cell Genomics."""

from ._version import get_versions

__version__ = get_versions()['version']
del get_versions
__maintainer__ = ', '.join([
    "Leander Dony",
    "David S. Fischer"
])
__author__ = ', '.join([
    "Leander Dony",
    "David S. Fischer"
])
__email__ = ', '.join([
    "leander.dony@helmholtz-muenchen.de",
    "david.fischer@helmholtz-muenchen.de"
])

import sfaira.data
import sfaira.genomes
import sfaira.models
import sfaira.train
import sfaira.interface as ui
