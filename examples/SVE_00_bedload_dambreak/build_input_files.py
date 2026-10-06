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

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=200,
        xa=0.,
        xb=10.,
        geo_file="MESH.geo",
        erodible_bed=True,
        zb=0.,
        zf=-0.1,
        verbose=False)

    # build initial condition:
    print("~~> create initial condition")
    build_ini_file(
        geo_file="MESH.geo", 
        ini_file="INIC.ini", 
        h_funct=h, 
        u=0., 
        verbose=False)
