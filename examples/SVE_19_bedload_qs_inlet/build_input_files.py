import os
import numpy as np
import matplotlib.pylab as plt
from arrakis import *

B = 0.2
Q = 0.0071

def zb(x):
    return -0.00427*x

def zf(x):
    return zb(x) - 1.5

def h(x):
    return 0.072

def u(x):
    return Q/(B*h(x))

if __name__ == '__main__':

    # build mesh:
    # ~~~~~~~~~~~
    print("~~> create mesh")
    build_geo_file(
        nx=150,
        xa=0.,
        xb=30.,
        geo_file="MESH_INI.geo",
        erodible_bed=False,
        variable_width=True,
        zb_funct=zb,
        zf_funct=zf,
        Bx=B,
        verbose=False)

    build_geo_file(
        nx=150,
        xa=0.,
        xb=30.,
        geo_file="MESH.geo",
        erodible_bed=True,
        variable_width=True,
        zb_funct=zb,
        zf_funct=zf,
        Bx=B,
        verbose=False)

    # build initial condition:
    # ~~~~~~~~~~~~~~~~~~~~~~~~
    print("~~> create initial condition")
    build_ini_file(
        geo_file="MESH.geo",
        ini_file="INIC.ini",
        h_funct=h,
        u_funct=u,
        verbose=False)

    # run init 
    # ~~~~~~~~
    print("~~> run init")
    os.system("arrakis.py bedload_qs_ini.yml")
    
    # build continuation file
    # ~~~~~~~~~~~~~~~~~~~~~~~
    print("~~> create initial condition (continuation)")
    os.system("cp RESU/RESfin.dat .")
    os.system("mv RESfin.dat CONT.ini")
    
    # Build boundary condition file:
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # water level function of time
    def q(t):
        return 0.0071

    def qs(t):
        qs_eq = 0.58e-05
        return qs_eq*(1. + 0.2*np.cos(2.*np.pi*t/1800.))

    # write .liq file
    print("~~> create boundary condition")
    build_bnd_file(
        bnd_file="BCL.liq", 
        ntimes=501, 
        t0=0., 
        tf=21600.,
        q_funct=q,
        qs_funct=qs,
        verbose=False)
        
    # plot
    bnd_data = np.loadtxt("BCL.liq", skiprows=3)
    fig, ax = plt.subplots(1, 1, figsize=(6.,4.5))
    ax.plot(bnd_data[:, 0], bnd_data[:, 2], color='k', ls='-')
    plt.show()
