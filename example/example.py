"""
Example script demonstrating how to call and run Xanthos.
"""

from xanthos.model import Xanthos


if __name__ == "__main__":

    # full path to parameterized config file
    ini = '/users/ladmin/repos/github/xanthos/example/config_hist.ini'

    # instantiate xanthos class
    xth = Xanthos(ini)

    # run intended configuration
    xth.execute()

    # clean up
    del xth