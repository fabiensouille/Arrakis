import numpy as np
import matplotlib.pylab as plt
from arrakis import *

def h(x):
    if x < 5.:
        h = 0.5
    else:
        h = 0.01
    return h

if __name__ == '__main__':

    nx = 10000

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=nx,
        xa=0.,
        xb=10.,
        geo_file="MESH_{}.geo".format(nx),
        erodible_bed=True,
        zb=0.,
        zf=-0.1,
        verbose=False)

    # build initial condition:
    print("~~> create initial condition")
    build_ini_file(
        geo_file="MESH_{}.geo".format(nx), 
        ini_file="INIC.ini", 
        h_funct=h, 
        u=0., 
        verbose=False)
