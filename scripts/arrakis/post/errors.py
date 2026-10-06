# -*- coding: utf-8 -*-
""" Error functions """
import numpy as np

def norm_l1(diff, mass=None):
    """
    Compute the l1 norm for the difference between v1 and v2

    @param diff (np.array) Vector containing (x1 - x2)
    @param mass (np.array) Mass factor (default: None)
    """
    if mass is None:
        return np.sum(abs(diff))/len(diff)
    else:
        return np.sum(mass*abs(diff))/np.sum(mass)

def norm_l2(diff, mass=None):
    """
    Compute the l2 norm for the difference between v1 and v2

    @param diff (np.array) Vector containing (x1 - x2)
    @param mass (np.array) Mass factor (default: None)
    """
    if mass is None:
        return np.sqrt(np.sum(abs(diff*diff))/len(diff))
    else:
        return np.sqrt(np.sum(mass*abs(diff*diff))/np.sum(mass))

def norm_linf(diff):
    """
    Compute the linf norm for the difference between v1 and v2

    @param diff (array) Vector containing (x1 - x2)
    """
    return np.max(abs(diff))
        
def compute_diff(data1, data2, relative=False):
    """
    Compute difference of two vectors

    @param data1 (array) usually the model or computation
    @param data2 (array) usually the reference data
    @param relative (bool) use relative difference
    """
    if relative:
        diff = data1 - data2
        for i, _ in enumerate(diff):
            if abs(diff[i]) > 1e-42:
                a_max = max(abs(data1[i]), abs(data2[i]))
                diff[i] = diff[i]/a_max
            else:
                diff[i] = 0.0
    else:
        diff = data1 - data2
    return diff
        
def compute_norm(v1, v2, norm='linf', mass=None, relative=False):
    """
    Compute norm of a vector

    @param diff (array) Vector containing (x1 - x2)
    @param mass (array) Mass factor (default: None)
    @param norm (str) type of norm
    """
    diff = compute_diff(v1, v2, relative=relative)
  
    if norm == 'l1':
        return norm_l1(diff, mass=mass)
    elif norm == 'l2':
        return norm_l2(diff, mass=mass)
    else:
        return norm_linf(diff)
