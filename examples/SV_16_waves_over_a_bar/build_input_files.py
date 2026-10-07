import numpy as np
import matplotlib.pylab as plt
from arrakis import *

PLOT = False

def zb(x):
    if 0.<=x and x<6.:
      return 0.
    elif 6.<=x and x<12.:
      return (x-6.)*(0.3/6.)
    elif 12.<=x and x<14.:
      return 0.3
    elif 14.<=x and x<17.:
      return 0.3-(x-14.)*(0.3/3.)
    elif 17.<=x and x<18.95:
      return 0.
    else:
      return (x-18.95)*(0.75/18.75)

def h(x):
     return max(0., 0.4 - zb(x))
      
if __name__ == '__main__':

    # Uniform mesh:
    # ~~~~~~~~~~~~ 
    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=500,
        xa=0.,
        xb=35.,
        geo_file="MESH.geo",
        zb_funct=zb,
        verbose=False)
        
    # plot
    if PLOT:
        mesh_data = np.loadtxt("MESH.geo", skiprows=2)
        fig, ax = plt.subplots(1, 1, figsize=(6.,4.5))
        ax.plot(mesh_data[:,0], mesh_data[:,1], color='k', ls='-', marker='s')
        plt.show()

    # build initial condition:
    print("~~> create initial condition")
    build_ini_file(
        geo_file="MESH.geo",
        ini_file="INIC.ini",
        h_funct=h,
        u=0.,
        verbose=False)

    # Build boundary condition file:
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # water level function of time
    def h(t):
        A = 0.029
        fm = 1.
        f1 = fm*(1. + 1./20.)
        f2 = fm*(1. - 1./20.)
        return 0.4 + A*(np.cos(2.*np.pi*f1*t) + np.cos(2.*np.pi*f2*t))

    # write .liq file
    print("~~> create boundary condition")
    build_bnd_file(
        bnd_file="BCL.liq", 
        ntimes=2001, 
        t0=0., 
        tf=100.,
        h_funct=h,
        verbose=False)

    # plot
    if PLOT:
        bnd_data = np.loadtxt("BCL.liq", skiprows=3)
        fig, ax = plt.subplots(1, 1, figsize=(6.,4.5))
        ax.plot(bnd_data[:, 0], bnd_data[:, 1], color='k', ls='-')
        plt.show()
