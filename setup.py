from setuptools import setup, find_packages
import versioneer

author = 'theislab'
author_email = 'david.fischer@helmholtz-muenchen.de'
description = "sfaira is a model and a data repository for single-cell data in a single python package. "

with open("README.rst", "r") as fh:
    long_description = fh.read()

setup(
    name='sfaira',
    author=author,
    author_email=author_email,
    description=description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(),
    install_requires=[
        'anndata>=0.7',
        'h5py',
        'numpy>=1.16.4',
        'pandas',
        'scipy>=1.2.1',
        'tqdm',
        'tensorflow>=2.0.0'  # TODO Remove and add to tensorflow profile
    ],
    extras_require={
        'tensorflow': [
            # TODO Add Tensorflow here again
        ],
        'kipoi': [
            'kipoi',
            'git-lfs'
        ],
        'plotting_deps': [
            "seaborn",
            "matplotlib",
            "sklearn"
        ],
        'data': [
            "scanpy",
            "loompy",
            "requests",
            "xlrd>=1.0.0"
        ],
        'extension': [
            "sfaira_extension",
        ],
        'docs': [
            'sphinx',
            'sphinx-autodoc-typehints',
            'sphinx_rtd_theme',
            'jinja2',
            'docutils',
        ],
    },
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
)
