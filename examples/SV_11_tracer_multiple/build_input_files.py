import numpy as np
import matplotlib.pylab as plt
from arrakis import *

def tracers(x, k):
    " tracers initial condition, k: index for tracer"
    # tracer 1
    if k==0:
        return np.exp(-2.5*(x-2.)**2)
    # tracer 2
    elif k==1:
        if x >= 4.5 and x <= 6.5:
            return 1.
        else:
            return 0.

if __name__ == '__main__':

    # build mesh:
    # ~~~~~~~~~~~
    print("~~> create mesh")
    build_geo_file(
        nx=100,
        xa=0.,
        xb=10.,
        geo_file="MESH.geo",
        verbose=False)

    # build initial condition:
    # ~~~~~~~~~~~~~~~~~~~~~~~
    print("~~> create initial condition")
    build_ini_file(
        geo_file="MESH.geo",
        ini_file="INIC.ini",
        h=1.,
        u=0.1,
        ntrac=2,
        t_funct=tracers,
        verbose=False)
        

    # Build boundary condition file:
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # water level function of time
    def h(t):
        return 1.

    def u(t):
        return 0.1

    def tracers_bc(t, k):
        if k==0:
            return 0.
        else:
            return 1.

    # write .liq file : case critical inlet
    print("~~> create boundary condition")
    build_bnd_file(
        bnd_file="BCL.liq",
        ntimes=11, 
        t0=0., 
        tf=100.,
        u_funct=u,
        ntrac=2,
        t_funct=tracers_bc,
        verbose=False)
