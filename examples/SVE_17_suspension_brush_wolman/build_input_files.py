import os
import numpy as np
import matplotlib.pylab as plt
from arrakis import *

"""
Lucien M Brush Jr and M Gordon Wolman. Knickpoint behavior in noncohesive material :, a
laboratory study. Geol. Soc. Am. Bull., 71(1) :59–74, 1960.
"""

PLOT = True

# parameters
I1 = -0.00125
I2 = -0.1
B = 0.21336
Q = 5.95e-4
q = Q/B
h = 0.02
u = q/h
print("u = ", u)

def zb(x):
    zb = 0.05 + I1*x
    if x>11. and x<11.3:
      zb = 0.05 + I1*11. + I2*(x-11.)
    elif x>11.3:
      zb = 0.05 + I1*11. + I2*(0.3) + I1*(x-11.3)
    return zb

def zf(x):
    return 0.05 + I1*x - 0.1
    
def h(x):
    return 0.05 + I1*x + 0.02 - zb(x)

if __name__ == '__main__':

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=250,
        xa=0.,
        xb=20.,
        geo_file="MESH_INI.geo",
        variable_width=True,
        erodible_bed=False,
        zb_funct=zb,
        Bx=B,
        verbose=False)

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=500,
        xa=0.,
        xb=20.,
        geo_file="MESH.geo",
        variable_width=True,
        erodible_bed=True,
        zb_funct=zb,
        zf_funct=zf,
        Bx=B,
        verbose=False)

    # build initial condition:
    print("~~> create initial condition")
    build_ini_file(
        geo_file="MESH.geo", 
        ini_file="INIC.ini", 
        h_funct=h,
        u=u,
        ntrac=1,
        t=0.024,
        verbose=False)
        
    if PLOT:
        # plot (x,z)
        mesh_data = np.loadtxt("MESH.geo", skiprows=2)
        ini_data = np.loadtxt("INIC.ini", skiprows=2, delimiter=',')
        x = mesh_data[:,0]
        zb = mesh_data[:,1] 
        zf = mesh_data[:,2] 
        h = ini_data[:,1]
        fig, ax = plt.subplots(1, 1, figsize=(5.,4.))
        ax.plot(x, zb+h, color='b', ls='-', lw=0.5, label='$h_0$')
        ax.plot(x, zb, color='k', ls='-', lw=0.5, label='$z_b$')
        ax.plot(x, zf, color='k', ls='--', lw=0.5, label='$z_f$')
        ax.fill_between(x, zb, zf, color='peru', alpha=0.5)
        ax.fill_between(x, np.min(zf), zf, color='saddlebrown', alpha=0.5)
        ax.fill_between(x, zb, zb+h, color='steelblue', alpha=0.25)
        plt.legend()
        ax.set_ylim([min(zf), zb[0]+h[0]])
        ax.set_xlim([0., 15.85])
        ax.set_ylabel("$z$ (m)")
        ax.set_xlabel("$x$ (m)")
        plt.savefig("FIG/geo_zb.png", dpi=300)
        plt.show()
        plt.close()

    # run init 
    print("~~> run init")
    os.system("arrakis.py brushwolman_ini.yml")
    
    # build continuation file
    print("~~> create initial condition (continuation)")
    os.system("cp RESU/RESfin.dat .")
    os.system("mv RESfin.dat CONT.ini")
