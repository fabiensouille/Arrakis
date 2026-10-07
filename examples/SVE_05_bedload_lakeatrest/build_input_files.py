import os
import numpy as np
import matplotlib.pylab as plt
from arrakis import *

def zb(x):
    if x <= 0.5:
      return np.sin(4.*np.pi*x)
    else:
      return np.sin(4.*np.pi*x) - 2.

def zf(x):
    return min(-1., zb(x))

def h(x):
    return max(1.2 - zb(x), 0.)

if __name__ == '__main__':

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=100,
        xa=0.,
        xb=1.,
        geo_file="MESH.geo",
        erodible_bed=True,
        zb_funct=zb,
        zf_funct=zf,
        verbose=False)

    # build initial condition:
    print("~~> create initial condition")
    build_ini_file(
        geo_file="MESH.geo", 
        ini_file="INIC.ini", 
        h_funct=h, 
        u=0.,
        verbose=False)
