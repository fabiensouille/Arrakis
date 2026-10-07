import numpy as np
import matplotlib.pylab as plt
from arrakis import *

PLOT = True

def zf(x):
    return 5. - (5./10000.)*x

def zb(x):
    x0 = 4000.
    x1 = 6000.
    if x>x0 and x<x1:
      return zf(x) + 1.
    else:
      return zf(x)

def h(x):
    return 5. + (zf(x)-zb(x))

if __name__ == '__main__':

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=250,
        xa=0.,
        xb=10000.,
        geo_file="MESH.geo",
        erodible_bed=True,
        zb_funct=zb,
        zf_funct=zf,
        verbose=False)

    # build initial condition:
    print("~~> create initial condition")

    build_ini_file(
        geo_file="MESH.geo",
        ini_file="INIC.ini",
        h_funct=h,
        q=10.,
        verbose=False)
        
    if PLOT:
        # plot (x,z)
        mesh_data = np.loadtxt("MESH.geo", skiprows=2)
        ini_data = np.loadtxt("INIC.ini", skiprows=2, delimiter=',')
        x = mesh_data[:,0]
        zb = mesh_data[:,1] 
        zf = mesh_data[:,2] 
        h = ini_data[:,1]
        print("cote retenue :", h[-1])
        fig, ax = plt.subplots(1, 1, figsize=(5.,4.))
        ax.plot(x, zb+h, color='b', ls='-', lw=0.5, label='$h_0$')
        ax.plot(x, zb, color='k', ls='-', lw=0.5, label='$z_b$')
        ax.plot(x, zf, color='k', ls='--', lw=0.5, label='$z_f$')
        ax.fill_between(x, zb, zf, color='peru', alpha=0.5)
        ax.fill_between(x, np.min(zf), zf, color='saddlebrown', alpha=0.5)
        ax.fill_between(x, zb, zb+h, color='steelblue', alpha=0.25)
        plt.legend()
        ax.set_ylim([min(zf), zb[0]+h[0]])
        ax.set_xlim([0., 10000.])
        ax.set_ylabel("$z$ (m)")
        ax.set_xlabel("$x$ (m)")
        plt.savefig("FIG/geo_zb.png", dpi=300)
        plt.show()
        plt.close()
