import os
import numpy as np
import matplotlib.pylab as plt
from arrakis import *

PLOT = False

def width(x):
    """ channel_width """
    width_entry = 0.064
    width_middle = 0.038
    pente = (width_entry-width_middle)/0.052
    if x < 0.3705:
        L = width_entry
    elif x >= 0.3705 and x <= 0.4225:
        L = width_entry - pente*(x-0.3705)
    elif x >= 0.4225 and x <= 0.4775:
        L = width_middle
    elif x >= 0.4775 and x <= 0.5295:
        L = width_middle + pente*(x-0.4775)
    else:
        L = width_entry
    return L

def zb(x):
    return 0.072

def zf(x):
    return 0.04

def h(x):
    return 0.066

def u(x):
    Q = 8.8333e-4
    return Q/(width(x)*h(x))

if __name__ == '__main__':

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=250,
        xa=0.,
        xb=0.9,
        geo_file="MESH_INI.geo",
        erodible_bed=False,
        variable_width=True,
        zb_funct=zb,
        zf_funct=zf,
        Bx_funct=width,
        verbose=False)

    build_geo_file(
        nx=250,
        xa=0.,
        xb=0.9,
        geo_file="MESH.geo",
        erodible_bed=True,
        variable_width=True,
        zb_funct=zb,
        zf_funct=zf,
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

    if PLOT = True:
        # plot (x,y)
        mesh_data = np.loadtxt("MESH.geo", skiprows=2)
        ini_data = np.loadtxt("INIC.ini", skiprows=2, delimiter=',')
        x = mesh_data[:,0]
        Lx = mesh_data[:,3]
        nx = len(x)
        fig, ax = plt.subplots(1, 1, figsize=(5.,4.))
        ax.plot(x, Lx/2., color='k', ls='-', lw=0.5,  label='$\pm L/2$')
        ax.plot(x,-Lx/2., color='k', ls='-', lw=0.5, label='')
        ax.plot(x, np.zeros(nx), color='k', lw=0.5, ls=':', label='')
        ax.fill_between(x, -Lx/2., Lx/2., color='steelblue', alpha=0.25)
        plt.legend()
        ax.set_ylim([-Lx[0]/2., Lx[0]/2.])
        ax.set_xlim([0., 0.9])
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
        ax.set_xlim([0., 0.9])
        ax.set_ylabel("$z$ (m)")
        ax.set_xlabel("$x$ (m)")
        plt.savefig("FIG/geo_zb.png", dpi=300)
        plt.show()
        plt.close()

    # run init 
    print("~~> run init")
    os.system("arrakis.py channel_constriction_ini.yml")
    
    # build continuation file
    print("~~> create initial condition (continuation)")
    os.system("cp RESU/RESfin.dat .")
    os.system("mv RESfin.dat CONT.ini")
