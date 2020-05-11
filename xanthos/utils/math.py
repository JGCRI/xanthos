"""
Mathematical helper functions.

Created on Feb 16, 2016
@author: lixi729

License:  BSD 2-Clause, see LICENSE and DISCLAIMER files

Copyright (c) 2017, Battelle Memorial Institute
"""

import logging
import numpy as np


def Size(l):
    """Get the size of a 2D array."""
    nrow = len(l)
    try:
        ncol = len(l[0])
    except:
        ncol = 1
    return nrow, ncol


def SizeR(l):
    nr = len(l)
    return nr


def SizeC(l):
    try:
        nc = len(l[0])
    except:
        nc = 1
    return nc


# Convert subscripts to linear indices
def sub2ind(arraySize, rowSub, colSub):
    """Convert subscripts to linear indices."""
    linearInd = []
    if len(rowSub) != len(colSub):
        logging.warning('def sub2ind at Rearranging: length of rowSub is not equal to length of colSub!')
    else:
        arr = tuple(arraySize)
        for i in range(0, len(rowSub)):
            temp = np.ravel_multi_index((rowSub[i], colSub[i]), arr, order='F')
            linearInd.append(temp)
    return np.array(linearInd)


# Convert linear indices to subscripts
def ind2sub(arraySize, index):
    """Convert linear indices to subscripts.

    :param index:   A list or 1d array
    """
    linearInd = np.zeros((len(index), 2), dtype=int)
    arr = tuple(arraySize)
    for i in range(0, len(index)):
        temp = np.unravel_index(index[i], arr, order='F')
        linearInd[i] = np.array(temp)
    return linearInd
