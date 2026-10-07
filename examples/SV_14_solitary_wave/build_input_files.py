import numpy as np
import matplotlib.pylab as plt
from arrakis import *

g = 9.81
H0 = 1.
l0 = np.sqrt(H0**2 + H0**3/0.35)

def h(x):
    return H0 + 0.35*(1./np.cosh((x-5.)/l0))**2

def u(x):
    c0 = l0*np.sqrt((g*h(x)**3)/(l0**2 - H0**2))
    return c0*(1. - 1./h(x))

if __name__ == '__main__':

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=500,
        xa=0.,
        xb=30.,
        geo_file="MESH.geo",
        zb=0.,
        verbose=False)

    # build initial condition:
    print("~~> create initial condition")
    build_ini_file(
        geo_file="MESH.geo",
        ini_file="INIC.ini",
        h_funct=h, 
        u_funct=u, 
        verbose=False)
