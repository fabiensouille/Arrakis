import os
import numpy as np
import matplotlib.pylab as plt
from arrakis import *

PLOT = False

def osc(x):
    return 5.*np.cos(1.-0.0002*x) + 0.5*np.cos(0.0035*x) + 0.2*np.sin(0.01*x)

def zb(x):
    return 30. - 0.0015*x - osc(x)

def zf(x):
    return 25. - 0.0015*x - 0.4*osc(x)

def h(x):
    return max(18.- zb(x), 28. - 0.0015*x - zb(x))

if __name__ == '__main__':

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=500,
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
        u=0.,
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

#    # build bnd condition:
#    h0 = 0.95
#    q0 = 0.1
#    a = 0.01
#    tm = 50.

#    def q_bnd(t):
#        if t < tm:
#            return q0 + a*t
#        else:
#            return q0 + 2.*a*tm -a*t

#    def h_bnd(t):
#        if t < tm:
#            return h0 + a*t
#        else:
#            return h0 + 2.*a*tm -a*t

#    print("~~> create boundary condition")
#    build_bnd_file(
#        bnd_file="BCQ.liq", 
#        ntimes=50, 
#        t0=0., 
#        tf=150.,
#        q_funct=q_bnd, 
#        verbose=False)
#        
#    build_bnd_file(
#        bnd_file="BCH.liq", 
#        ntimes=50, 
#        t0=0., 
#        tf=150.,
#        h_funct=h_bnd, 
#        verbose=False)
        
    # run init 
    print("~~> run init")
    os.system("arrakis.py bedload_stream_ini.yml")
    
    # build continuation file
    print("~~> create initial condition (continuation)")
    os.system("cp RESU/RESfin.dat .")
    os.system("mv RESfin.dat CONT.ini")
