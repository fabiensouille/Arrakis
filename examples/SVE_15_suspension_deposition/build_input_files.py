import os
import numpy as np
import matplotlib.pylab as plt
from arrakis import *

PLOT = False

if __name__ == '__main__':

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=100,
        xa=0.,
        xb=10.,
        geo_file="MESH.geo",
        erodible_bed=True,
        zb=0.0,
        zf=-0.5,
        verbose=False)

    # build initial condition:
    print("~~> create initial condition")
    build_ini_file(
        geo_file="MESH.geo", 
        ini_file="INIC.ini", 
        h=1.,
        u=0.,
        ntrac=1,
        t=1.,
        verbose=False)
