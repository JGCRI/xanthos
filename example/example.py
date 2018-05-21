"""
Example script demonstrating how to call and run Xanthos.
"""

from xanthos.model import Xanthos


def run(ini):

    # instantiate model
    xth = Xanthos(ini)

    # run intended configuration
    xth.execute()

    return xth


if __name__ == "__main__":

    # full path to parameterized config file
    ini = '/Users/d3y010/repos/github/xanthos/example/pm_abcd_mrtm.ini'

    # run the model
    xth = run(ini)

    del xth