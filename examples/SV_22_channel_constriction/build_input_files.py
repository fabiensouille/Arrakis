import numpy as np
import matplotlib.pylab as plt
from arrakis import *
"""
Garegnani, G.; Rosatti, G. & Bonaventura, L. On the range of validity of the Exner-based models for mobile-bed river flow simulations Journal of Hydraulic Research, Taylor & Francis, 2013, 51, 380-391
"""
PLOT = False

def width(x):
    """ channel_width """
    width_entry = 14.
    width_middle = 6.
    pente = (width_entry-width_middle)/200.
    if x < 500.:
        L = 14.
    elif x >= 500. and x <= 700.:
        L = width_entry - pente*(x-500.)
    elif x >= 700. and x <= 900.:
        L = width_middle
    elif x >= 900. and x <= 1100.:
        L = width_middle + pente*(x-900.)
    else:
        L = 14.
    return L

def zb(x):
    pente = 0.5/100.
    zb = 5. -pente*(x-2000.)
    return zb

def h(x):
    return 10.
    
def u(x):
    Q = 400.
    return Q/(width(x)*h(x))

if __name__ == '__main__':

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=200,
        xa=0.,
        xb=2000.,
        geo_file="MESH.geo",
        variable_width=True,
        zb_funct=zb,
        Bx_funct=width,
        verbose=False)

    # build initial condition:
    print("~~> create initial condition")
    build_ini_file(
        geo_file="MESH.geo",
        ini_file="INIC.ini",
        h_funct=h, 
        u_funct=u, 
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
        ax.set_xlim([0., 2000.])
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
        ax.set_xlim([0., 2000.])
        ax.set_ylabel("$z$ (m)")
        ax.set_xlabel("$x$ (m)")
        plt.savefig("FIG/geo_zb.png", dpi=300)
        plt.show()
        plt.close()
