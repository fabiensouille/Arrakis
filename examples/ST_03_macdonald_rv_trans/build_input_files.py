import numpy as np
import matplotlib.pylab as plt
from scipy.integrate import quad
from scipy.misc import derivative
from arrakis import *
"""
03. Mac Donald solution: 
************************ 
Description: subcritical flow in a rectangular channel 
             with variable width
example : Channel Contraction Test Problem of the article:
Xue, Z.; Zhou, L. & Liu, D. Accurate Numerical Modeling for 1D 
Open-Channel Flow with Varying Topography Water, MDPI, 2023
see also: MacDonald, I.; Baines, M.J.; Nichols, N.K.; Samuels, P.G.
Analytic benchmark solution of open-channel flows. J.Hydraul.Eng. 1997
"""
# parameters
g = 9.80665 # Gravity
Q = 20.     # Flow rate
n = 0.03    # Manning coefficient

PLOT = False

def width(x):
    """ channel_width """
    xa = x/1000.
    L = 10. - 64.*( xa**2 -2.*xa**3 + xa**4 )
    return L

def depth(x):
    """ channel_depth """
    xa = x/1000.
    if x<500.:
       H = -(1./40.) + (1./(1. + 2.*(xa-0.5)**2))
    else:
        a0 = 1.5
        ak = [-0.230680, 0.248267, -0.228271]
        H = a0*np.exp((x/4000.)-(1./4.))
        for k in range(3):
          H += ak[k]*np.exp(15.*(k+1) - 30.*(k+1)*xa)
    return H

def slope(x):
    """ channel_slope """
    L = width(x)
    H = depth(x)
    Q2 = Q**2
    L2 = L**2
    L3 = L**3
    H2 = H**2
    H3 = H**3
    dL = derivative(width, x, dx=0.00001)
    dH = derivative(depth, x, dx=0.00001)
    u2 = Q2/(H2*L2)
    aux = ((2./L) + (1./H))
    Sf = (u2*n**2)*( np.sign(aux)*(np.abs(aux))**(4/3) )
    S0 = Sf + (1. - u2/(g*H))*dH - (u2/g)*(dL/L)
    return S0

def bottom(x, xb=1000.):
    """ channel_bottom """
    zb = quad(func=slope, a=x, b=xb)
    zb = zb[0]
    return zb

def velocity(x):
    return Q/max(0., width(x)*depth(x))

if __name__ == '__main__':

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=150,
        xa=0.,
        xb=1000.,
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
        ax.set_xlim([0., 1000.])
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
        ax.set_xlim([0., 1000.])
        ax.set_ylabel("$z$ (m)")
        ax.set_xlabel("$x$ (m)")
        plt.savefig("FIG/geo_zb.png", dpi=300)
        plt.show()
        plt.close()

