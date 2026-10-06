import numpy as np
import matplotlib.pylab as plt
from arrakis import *

def h(x):
    if x > 5.:
        h = 0.2
    else:
        h = 1.
    return h

if __name__ == '__main__':

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=100,
        xa=0.,
        xb=10.,
        geo_file="MESH.geo", 
        verbose=False)

    # build initial condition:
    print("~~> create initial condition")
    build_ini_file(
        geo_file="MESH.geo", 
        ini_file="INIC.ini", 
        h_funct=h, 
        u=0., 
        verbose=False)
