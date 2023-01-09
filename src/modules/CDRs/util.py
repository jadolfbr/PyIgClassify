import math


def get_norm_distance(length, distance):
    """
    Calculate the normalized distance

    :param length: int
    :param distance: float
    :rtype: float
    """
    return distance/(length*2)

def get_norm_distance_deg(norm_distance):
    """
    Get the normalized distance in degrees.  This is probably what you want. If it is some crazy number, will return 360.
    :param norm_distance:
    :rtype: float
    """
    if norm_distance == -1.0:
        return -1.0
    elif int(norm_distance) >= 25:
        return 360
    elif norm_distance:
        return math.acos(1.0 - ( float( norm_distance )/2.0)) * (180/ math.pi)
    else:
        return 0