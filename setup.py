"""
@author   Chris R. Vernon
@email:   chris.vernon@pnnl.gov
@Project: Xanthos 2.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""
import sys


class VersionError(Exception):
    def __init__(self, *args, **kwargs):
        Exception.__init__(self, *args, **kwargs)


py_version = '{0}.{1}'.format(sys.version_info.major, sys.version_info.minor)
if py_version != '2.7':
    raise VersionError("Xanthos must be ran using Python 2.7.  You are using Python {0}.".format(py_version))


try:
    from setuptools import setup, find_packages
except ImportError:
    raise("Must have setuptools installed to run setup.py. Please install and try again.")


def readme():
    with open('README.md') as f:
        return f.read()


def get_requirements():
    with open('requirements.txt') as f:
        return f.read().split()


setup(
    name='xanthos',
    version='2.0.0',
    packages=find_packages(),
    url='https://github.com/jgcri/xanthos',
    license='BSD 2-Clause',
    author='Chris R. Vernon; Xinya Li',
    author_email='chris.vernon@pnnl.gov; xinya.li@pnl.gov',
    description='A global hydrologic model for GCAM',
    long_description=readme(),
    install_requires=get_requirements()
)
