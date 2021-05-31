import os

from setuptools import setup, find_packages
import versioneer

author = 'theislab'
author_email = 'david.fischer@helmholtz-muenchen.de'
description = "sfaira is a model and a data repository for single-cell data in a single python package."

with open("README.rst", "r") as fh:
    long_description = fh.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()


def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths


WD = os.path.dirname(__file__)
templates = package_files(os.path.join(WD, "sfaira", "commands", "templates"))

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
    packages=find_packages(include=['sfaira', 'sfaira.*']),
    package_data={'': templates},
    entry_points={
        'console_scripts': [
            'sfaira=sfaira.cli:main',
        ],
    },
    install_requires=requirements,
    extras_require={
        'tensorflow': [
            # 'tensorflow>=2.0.0'  # TODO Add Tensorflow here again
        ],
        'plotting_deps': [
            "seaborn",
            "matplotlib",
            "sklearn"
        ],
        'extension': [
            "sfaira_extension",
        ],
    },
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
)
