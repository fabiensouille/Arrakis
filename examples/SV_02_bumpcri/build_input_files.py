import numpy as np
import matplotlib.pylab as plt
from arrakis import *

def zb(x):
    return 0.1 + 0.1*np.exp(-(x-5.)**2.)

def h(x):
    return max(0.4 - zb(x), 0.)

if __name__ == '__main__':

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=100,
        xa=0.,
        xb=10.,
        geo_file="MESH.geo",
        zb_funct=zb,
        verbose=False)

    # build initial condition:
    print("~~> create initial condition")
    build_ini_file(
        geo_file="MESH.geo",
        ini_file="INIC.ini",
        h_funct=h, 
        q=0.15, 
        verbose=False)
