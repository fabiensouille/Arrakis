import numpy as np
import matplotlib.pylab as plt
from arrakis import *

def zb(x):
    if 1.4 < x and x < 1.6:
      return 2. + 0.25*(np.cos(10.*np.pi*(x-0.5)) + 1.)
    else:
      return 2.

def h(x):
    if 1.1 < x and x < 1.2:
      return max(0., 3. - zb(x) + 0.001)
    else:
      return max(0., 3. - zb(x))

if __name__ == '__main__':

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=80,
        xa=0.,
        xb=2.,
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
