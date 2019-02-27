"""
@author   Chris R. Vernon
@email:   chris.vernon@pnnl.gov
@Project: Xanthos 2.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

from xanthos.install_supplement import InstallSupplement


class VersionError(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)


try:
    from setuptools import setup, find_packages
except ImportError:
    print("Must have setuptools installed to run setup.py. Please install and try again.")
    raise


def readme():
    with open('README.md') as f:
        return f.read()


def get_requirements():
    with open('requirements.txt') as f:
        return f.read().split()


# install supplemental and example data from Zenodo
InstallSupplement().fetch_zenodo()


setup(
    name='xanthos',
    version='2.3.0',
    packages=find_packages(),
    url='https://github.com/jgcri/xanthos',
    license='BSD 2-Clause',
    author='Chris R. Vernon; Xinya Li',
    author_email='chris.vernon@pnnl.gov; xinya.li@pnl.gov',
    description='A global hydrologic model for GCAM',
    long_description=readme(),
    install_requires=get_requirements(),
    python_requires='>=2.7, !=3.0.*, !=3.1.*, !=3.2.*, <4',
)
