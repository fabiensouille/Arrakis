import numpy as np
import matplotlib.pylab as plt
from arrakis import *

PLOT = False

def zb(x):
    if 0.<= x and x < 10.:
        return 0.
    else:
        return (0.5/10.)*(x-10.)

def h(x):
     return max(0., 0.47 - zb(x))
      
if __name__ == '__main__':

    # Uniform mesh for run-up case :
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=500,
        xa=0.,
        xb=25.,
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

    # Non-uniform mesh for run-up2 case:
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    def mesh_density(x):
        if 0.<=x and x<2.5:
            return 0.5
        elif 5.<=x and x<15.:
            return 2.
        else:
            return 1.

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=500,
        xa=0.,
        xb=25.,
        geo_file="MESH_2.geo",
        density=mesh_density,
        zb_funct=zb,
        verbose=False)

    # plot
    if PLOT:
        mesh_data = np.loadtxt("MESH_2.geo", skiprows=2)
        fig, ax = plt.subplots(1, 1, figsize=(6.,4.5))
        ax.plot(mesh_data[:,0], mesh_data[:,1], color='k', ls='-', marker='s')
        plt.show()

    # build initial condition:
    print("~~> create initial condition")
    build_ini_file(
        geo_file="MESH_2.geo",
        ini_file="INIC_2.ini",
        h_funct=h, 
        u=0., 
        verbose=False)

    # Build boundary condition file:
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # water level function of time
    def zs(t):
        A = 0.015
        fm = 1.
        f1 = fm*(1. + 1./20.)
        f2 = fm*(1. - 1./20.)
        return 0.47 + A*(np.cos(2.*np.pi*f1*t) + np.cos(2.*np.pi*f2*t))

    # write .liq file
    print("~~> create boundary condition")
    build_bnd_file(
        bnd_file="BCL.liq", 
        ntimes=2001, 
        t0=0., 
        tf=100.,
        zs_funct=zs, 
        verbose=False)

    # plot
    if PLOT:
        bnd_data = np.loadtxt("BCL.liq", skiprows=3)
        fig, ax = plt.subplots(1, 1, figsize=(6.,4.5))
        ax.plot(bnd_data[:, 0], bnd_data[:, 1], color='k', ls='-')
        plt.show()
