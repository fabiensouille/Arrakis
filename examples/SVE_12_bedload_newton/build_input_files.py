import os
import numpy as np
import matplotlib.pylab as plt
from arrakis import *

B = 0.3048
Q = 0.00566

def zb(x):
    return -0.00416*x

def zf(x):
    return zb(x) - 1.5

def h(x):
    return 0.041

def u(x):
    return Q/(B*h(x))

if __name__ == '__main__':

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=100,
        xa=0.,
        xb=9.14,
        geo_file="MESH_INI.geo",
        erodible_bed=False,
        variable_width=True,
        zb_funct=zb,
        zf_funct=zf,
        Bx=B,
        verbose=False)

    build_geo_file(
        nx=100,
        xa=0.,
        xb=9.14,
        geo_file="MESH.geo",
        erodible_bed=True,
        variable_width=True,
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
        u_funct=u,
        verbose=False)

    # run init 
    print("~~> run init")
    os.system("arrakis.py bedload_newton_ini.yml")
    
    # build continuation file
    print("~~> create initial condition (continuation)")
    os.system("cp RESU/RESfin.dat .")
    os.system("mv RESfin.dat CONT.ini")
    
    # Bedload bnd condition
    q = Q/B
    rhos = 2650.
    Cs = 0.88
    qs = q*Cs/rhos
    print("~~> Bedload bnd condition ")
    print("qs = ", qs)
