from setuptools import setup, find_packages
import versioneer

author = 'theislab'
author_email = 'david.fischer@helmholtz-muenchen.de'
description = "sfaira is a model and a data repository for single-cell data in a single python package."

with open("README.rst", "r") as fh:
    long_description = fh.read()

setup(
    name='sfaira',
    author=author,
    author_email=author_email,
    description=description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
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
            # 'tensorflow>=2.0.0'  # TODO Add Tensorflow here again
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
            "xlrd==1.*",
            "openpyxl",
        ],
        'extension': [
            "sfaira_extension",
        ],
    },
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
)
