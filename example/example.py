"""
Example script demonstrating how to call and run Xanthos.
"""
import os

from xanthos import Xanthos


def run(ini):

    # instantiate model
    xth = Xanthos(ini)

    # run intended configuration
    xth.execute()

    return xth


if __name__ == "__main__":

    # full path to parameterized config file
    ini = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'pm_abcd_mrtm.ini')

    # run the model
    xth = run(ini)

    del xth