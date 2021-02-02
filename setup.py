import importlib
import os

from setuptools import setup, find_packages
import versioneer

author = 'theislab'
author_email = 'david.fischer@helmholtz-muenchen.de'
description = "sfaira is a model and a data repository for single-cell data in a single python package."
sfaira_module_path_init_py = importlib.util.find_spec("sfaira")
sfaira_module_path = sfaira_module_path_init_py.origin[:-12]


def walker(base, *paths):
    file_list = set([])
    cur_dir = os.path.abspath(os.curdir)

    os.chdir(base)
    try:
        for path in paths:
            for dname, dirs, files in os.walk(path):
                for f in files:
                    file_list.add(os.path.join(dname, f))
    finally:
        os.chdir(cur_dir)

    return list(file_list)


with open("README.rst", "r") as fh:
    long_description = fh.read()

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

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
    package_data={
        'sfaira': walker(
            os.path.dirname(sfaira_module_path),
            'create/templates'  # TODO LH
        ),
    },
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
