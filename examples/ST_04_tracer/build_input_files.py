import numpy as np
import matplotlib.pylab as plt
from arrakis import *

def tracers(x, k):
    " tracers initial condition, k: index for tracer"
    if x >= 4.5 and x <= 6.5:
        d = 1.
    else:
        d = 0.
    return d + np.exp(-2.5*(x-2.)**2)

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
        h=1.,
        u=0.1,
        ntrac=1,
        t_funct=tracers,
        verbose=False)
