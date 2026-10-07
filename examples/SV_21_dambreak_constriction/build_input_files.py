import numpy as np
import matplotlib.pylab as plt
from arrakis import *
"""
----------------------------------------------------------------
05. Channel with constriction (dambreak) 
*****************************
Description: dambreak with experimental measurments
section 6.4. of : A kinetic interpretation of the section-averaged 
Saint-Venant system for natural river hydraulics
see also: Goutal, N. & Maurel, F. A finite volume solver for 1D 
shallow-water equations applied to an actual river International 
Journal for Numerical Methods in Fluids, 2002
"""
def width(x):
    """ channel_width """
    width_entry = 0.5
    width_middle = 0.1
    pente = (width_entry-width_middle)/0.2
    if x < 14.:
        L = 0.5
    elif x >= 14. and x <= 14.2:
        L = width_entry - pente*(x-14.)
    elif x >= 14.2 and x <= 15.2:
        L = width_middle
    elif x >= 15.2 and x <= 15.4:
        L = width_middle + pente*(x-15.2)
    else:
        L = 0.5
    return L

def h(x):
    if x < 6.1:
      h = 0.3
    else:
      h = 0.003
    return h

if __name__ == '__main__':

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=400,
        xa=0.,
        xb=19.3,
        geo_file="MESH.geo",
        variable_width=True,
        zb=0.,
        Bx_funct=width,
        verbose=False)

    # build initial condition:
    print("~~> create initial condition")
    build_ini_file(
        geo_file="MESH.geo",
        ini_file="INIC.ini",
        h_funct=h, 
        u=0., 
        verbose=False)
