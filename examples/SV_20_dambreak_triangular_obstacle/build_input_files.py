import numpy as np
import matplotlib.pylab as plt
from arrakis import *
"""
04. Channel with constriction (dambreak) 
*****************************
Description: dambreak with experimental measurments
section 6.4. of : A kinetic interpretation of the section-averaged 
Saint-Venant system for natural river hydraulics
see also: Goutal, N. & Maurel, F. A finite volume solver for 1D 
shallow-water equations applied to an actual river International 
Journal for Numerical Methods in Fluids, 2002
"""
def zb(x):
    pente = 0.4/3
    if x < 25.5:
      zb = 0.
    elif x >= 25.5 and x <= 28.5:
      zb = pente*(x-25.5)
    elif x >= 28.5 and x <= 31.5:
      zb = 0.4 - pente*(x-28.5)
    else:
      zb = 0.
    return zb

def h(x):
    if x < 15.5:
      h = 0.75
    else:
      h = 0.
    return h

if __name__ == '__main__':

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=250,
        xa=0.,
        xb=38.,
        geo_file="MESH.geo",
        variable_width=True,
        zb_funct=zb,
        Bx=1.75,
        verbose=False)

    # build initial condition:
    print("~~> create initial condition")
    build_ini_file(
        geo_file="MESH.geo",
        ini_file="INIC.ini",
        h_funct=h, 
        u=0., 
        verbose=False)
