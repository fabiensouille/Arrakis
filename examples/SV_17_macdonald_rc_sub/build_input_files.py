import numpy as np
import matplotlib.pylab as plt
from scipy.integrate import quad
from arrakis import *
"""
01. Mac Donald solution: 
************************
Description: subcritical flow in a rectangular channel 
             with constant width
example 1 of : A kinetic interpretation of the section-averaged 
Saint-Venant system for natural river hydraulics
"""
# parameters
g = 9.80665 # Gravity
Q = 20.     # Flow rate
B = 10.     # width
n = 0.03    # Manning coefficient

PLOT = False

def depth(x):
    fact = -16.*((x/1000.) - (1./2.))**2
    H = ((4./g)**(1./3.))*( 1. + 0.5*np.exp(fact) )
    return H

def depth_derivative(x):
    fact = -16.*((x/1000.) - (1./2.))**2
    dfact = -(32./1000.)*((x/1000.) - (1./2.))
    dHdx = ((4./g)**(1./3.))*0.5*dfact*np.exp(fact)
    return dHdx

def slope(x):
    H = depth(x)
    dH = depth_derivative(x)
    dz = (1. - (4./(g*H**3)))*dH \
       + 0.36*( ((10. + 2*H)**(4./3.))/((10.*H)**(10./3.)) )
    return dz

def bottom(x, xb=1000.):
    zb = quad(func=slope, a=x, b=xb)
    zb = zb[0]
    return zb

def velocity(x):
    return Q/(B*depth(x))

if __name__ == '__main__':

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=150,
        xa=0.,
        xb=1000.,
        geo_file="MESH.geo",
        zb_funct=bottom,
        variable_width=True,
        Bx=B,
        verbose=False)
        
    # build initial condition:
    print("~~> create initial condition")
    build_ini_file(
        geo_file="MESH.geo",
        ini_file="INIC.ini",
        h_funct=depth,
        u_funct=velocity,
        verbose=False)

    # plot
    if PLOT:
        mesh_data = np.loadtxt("MESH.geo", skiprows=2)
        ini_data = np.loadtxt("INIC.ini", skiprows=2)
        x = mesh_data[:,0]
        zb = mesh_data[:,1] 
        h = ini_data[:,1]
        fig, ax = plt.subplots(1, 1, figsize=(5.,4.))
        ax.plot(x, zb+h, color='b', ls='-', lw=0.5, label='$h_0$')
        ax.plot(x, zb, color='k', ls='-', lw=0.5, label='$z_b$')
        ax.fill_between(x, np.min(zb), zb, color='grey', alpha=0.5)
        ax.fill_between(x, zb, zb+h, color='steelblue', alpha=0.25)
        plt.legend()
        ax.set_ylim([zb[-1], zb[0]+h[0]])
        ax.set_xlim([0., 1000.])
        ax.set_ylabel("$z$ (m)")
        ax.set_xlabel("$x$ (m)")
        plt.savefig("FIG/geo_zb.png", dpi=300)
        plt.show()
        plt.close()

