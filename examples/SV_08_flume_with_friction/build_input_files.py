import numpy as np
import matplotlib.pylab as plt
from arrakis import *

def zb(x):
    return 5. - (5./10000.)*x

if __name__ == '__main__':

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=100,
        xa=0.,
        xb=10000.,
        geo_file="MESH.geo",
        zb_funct=zb,
        verbose=False)

    # build initial condition:
    print("~~> create initial condition")

    build_ini_file(
        geo_file="MESH.geo",
        ini_file="INIC.ini",
        h=5.,
        q=10.,
        verbose=False)
