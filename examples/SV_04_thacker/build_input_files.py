import numpy as np
import matplotlib.pylab as plt
from arrakis import *

def zb(x):
    return 0.5*((x-5.)**2/3.**2-1.)

def h(x):
    return max(0., -0.5*(((x-5.)/3.+1./6.)**2-1.))

if __name__ == '__main__':

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=200,
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
        u=0., 
        verbose=False)
