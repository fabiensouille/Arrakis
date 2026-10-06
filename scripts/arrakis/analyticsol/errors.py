import numpy as np

def error_L2(array, array_ref, absolute=True, eps=1.e-16):
    """
    Computes error in L2 norm between an array and a reference
    """
    error = 0
    npoin = len(array)
    for i in range(npoin):
        if absolute:
            error += (array_ref[i] - array[i])**2
        else:
            error += ((array_ref[i] - array[i])/(eps+array_ref[i]))**2
    error = np.sqrt(error/npoin)
    return error

def error_L1(array, array_ref, absolute=True, eps=1.e-16):
    """
    Computes error in L1 norm between an array and a reference
    """
    error = 0
    npoin = len(array)
    for i in range(npoin):
        if absolute:
            error += abs(array_ref[i] - array[i])
        else:
            error += abs((array_ref[i] - array[i]))/(eps+array_ref[i])
    error = error/npoin
    return error
