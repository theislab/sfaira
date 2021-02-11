# -*- coding: utf-8 -*-
"""A Data and Model Zoo for Single-Cell Genomics."""

import warnings

# import sfaira.interface as ui
# import sfaira.train
# import sfaira.models
import sfaira.genomes
import sfaira.data
import sfaira.consts as consts
import sfaira.versions
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


# Monkey-patching warnings to prettify warning messages
def custom_formatwarning(msg, *args, **kwargs):
    # ignore everything except the message
    return str(msg) + '\n'


warnings.formatwarning = custom_formatwarning
