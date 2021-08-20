"""
Paths to cache directories used throughout the code.
"""

import os

CACHE_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(__file__))), "cache")

CACHE_DIR_DATABASES = os.path.join(CACHE_DIR, "dataset_meta")
CACHE_DIR_DATABASES_CELLXGENE = os.path.join(CACHE_DIR_DATABASES, "cellxgene")

CACHE_DIR_GENOMES = os.path.join(CACHE_DIR, "genomes")

CACHE_DIR_ONTOLOGIES = os.path.join(CACHE_DIR, "ontologies")
