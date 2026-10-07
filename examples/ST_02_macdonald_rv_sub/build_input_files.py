import numpy as np
import matplotlib.pylab as plt
from scipy.integrate import quad
from arrakis import *
"""
02. Mac Donald solution: 
************************ 
Description: subcritical flow in a rectangular channel 
             with variable width
example 2 of : A kinetic interpretation of the section-averaged 
Saint-Venant system for natural river hydraulics
see also: Mac Donald 1995. Comparison of some steady state 
Saint-Venant solvers for some test problems with analytic solutions
"""
# parameters
g = 9.80665 # Gravity
Q = 20.     # Flow rate
B = 10.     # width
n = 0.03    # Manning coefficient

PLOT = False

def width(x):
    """ channel_width """
    fact1 = -64.*((x/5000.) - (1./3.))**2
    fact2 = -64.*((x/5000.) - (2./3.))**2
    L = 10. - 5.*np.exp(fact1) - 5.*np.exp(fact2)
    return L

def depth(x):
    """ channel_depth """
    fact1 = -64.*((x/5000.) - (1./3.))**2
    fact2 = -64.*((x/5000.) - (2./3.))**2
    H = 1. + 0.5*np.exp(fact1) + 0.5*np.exp(fact2)
    return H

def width_derivative(x):
    """ channel_width """
    fact1 = -64.*((x/5000.) - (1./3.))**2
    fact2 = -64.*((x/5000.) - (2./3.))**2
    dfact1 = -(128./5000.)*((x/5000.) - (1./3.))
    dfact2 = -(128./5000.)*((x/5000.) - (2./3.))
    dLdx = - 5.*dfact1*np.exp(fact1) - 5.*dfact2*np.exp(fact2)
    return dLdx

def depth_derivative(x):
    """ channel_depth """
    fact1 = -64.*((x/5000.) - (1./3.))**2
    fact2 = -64.*((x/5000.) - (2./3.))**2
    dfact1 = -(128./5000.)*((x/5000.) - (1./3.))
    dfact2 = -(128./5000.)*((x/5000.) - (2./3.))
    dHdx = - 0.5*dfact1*np.exp(fact1) - 0.5*dfact2*np.exp(fact2)
    return dHdx

def slope(x):
    """ channel_slope """
    L = width(x)
    H = depth(x)
    dL = width_derivative(x)
    dH = depth_derivative(x)
    dz = (1. - (Q**2)/(g*(L**2)*(H**3)))*dH \
       + (Q**2)*(n**2)*( ((2*H + L)**(4./3.))/((L*H)**(10./3.)) ) \
       - (Q**2)*( (dL)/(g*(L**3)*(H**2)) )
    return dz

def bottom(x, xb=5000.):
    """ channel_bottom """
    zb = quad(func=slope, a=x, b=xb)
    zb = zb[0]
    return zb
    
def velocity(x):
    return Q/(width(x)*depth(x))

if __name__ == '__main__':

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=150,
        xa=0.,
        xb=5000.,
        geo_file="MESH.geo",
        variable_width=True,
        zb_funct=bottom,
        Bx_funct=width,
        verbose=False)
        
    # build initial condition:
    print("~~> create initial condition")
    build_ini_file(
        geo_file="MESH.geo",
        ini_file="INIC.ini",
        h_funct=depth,
        u_funct=velocity,
        verbose=False)

    if PLOT:
        # plot (x,y)
        mesh_data = np.loadtxt("MESH.geo", skiprows=2)
        ini_data = np.loadtxt("INIC.ini", skiprows=2)
        x = mesh_data[:,0]
        Lx = mesh_data[:,2]
        nx = len(x)
        fig, ax = plt.subplots(1, 1, figsize=(5.,4.))
        ax.plot(x, Lx/2., color='k', ls='-', lw=0.5,  label='$\pm L/2$')
        ax.plot(x,-Lx/2., color='k', ls='-', lw=0.5, label='')
        ax.plot(x, np.zeros(nx), color='k', lw=0.5, ls=':', label='')
        ax.fill_between(x, -Lx/2., Lx/2., color='steelblue', alpha=0.25)
        plt.legend()
        ax.set_ylim([-Lx[0]/2., Lx[0]/2.])
        ax.set_xlim([0., 5000.])
        ax.set_ylabel("$y$ (m)")
        ax.set_xlabel("$x$ (m)")
        plt.savefig("FIG/geo_width.png", dpi=300)
        plt.show()
        plt.close()

        # plot (x,z)
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
        ax.set_xlim([0., 5000.])
        ax.set_ylabel("$z$ (m)")
        ax.set_xlabel("$x$ (m)")
        plt.savefig("FIG/geo_zb.png", dpi=300)
        plt.show()
        plt.close()

