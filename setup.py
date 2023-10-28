"""
@author   Chris R. Vernon
@email:   chris.vernon@pnnl.gov
@Project: Xanthos 2.0

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

from setuptools import setup, find_packages


def readme():
    with open('README.md') as f:
        return f.read()


setup(
    name='xanthos',
    version='2.4.1',
    packages=find_packages(),
    url='https://github.com/JGCRI/xanthos',
    license='BSD2-Simplified',
    author='Chris R. Vernon',
    author_email='chris.vernon@pnnl.gov',
    description='A global hydrologic modeling framework',
    long_description=readme(),
    long_description_content_type="text/markdown",
    install_requires=[
        "numpy>=1.24.4",
        "scipy>=1.6",
        "pandas~=2.0.3",
        "configobj>=5.0.8",
        "joblib>=1.3.2",
        "matplotlib>=3.7.2",
        "xarray>=2023.8.0",
        "requests"
    ],
    python_requires='>=3.6',
    include_package_data=True
)
