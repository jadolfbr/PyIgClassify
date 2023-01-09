import math
import numpy

def distance_numpy(array1, array2):
    """
    Get the distance between two points
    :param array1: numpy.Array
    :param array2: numpy.Array
    :rtype: float
    """

    return numpy.linalg.norm(array1-array2)


def distance(x1, y1, z1, x2, y2, z2):
    """
    Get the distance between variables.
    :param x1: float
    :param y1: float
    :param z1: float
    :param x2: float
    :param y2: float
    :param z2: float
    :rtype: float
    """

    return math.sqrt( math.pow( (x1-x2), 2) ) + math.sqrt( math.pow( (y1-y2), 2) )  + math.sqrt( math.pow( (z1-z2), 2) )