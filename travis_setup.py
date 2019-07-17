#!/usr/bin/env python

import pkg_resources
from xanthos.install_supplement import InstallSupplement


def install_supplement():
    """Install supplement for use in Travis-CI testing."""

    # install example supplement for testing in Travis-CI
    InstallSupplement(pkg_resources.resource_filename('xanthos', 'test'))


if __name__ == '__main__':

    install_supplement()
