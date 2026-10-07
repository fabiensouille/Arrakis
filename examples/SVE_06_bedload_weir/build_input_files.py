import os
import numpy as np
import matplotlib.pylab as plt
from arrakis import *

xa = 8.
xb = 9.
xc = 10.
xd = 11.
zbot = 0.
ztop = 1.

def zb(x):
    pente = (xb-xa/ztop-zbot)
    if x<xa:
        zb = zbot
    elif xa<=x and x<xb:
        zb = zbot + pente*(x-xa)
    elif xb<=x and x<xc:
        zb = ztop
    elif xc<=x and x<xd:
        zb = ztop - pente*(x-xc)
    else:
        zb = zbot
    if x<0.5*(xb+xc):
        zb = max(zb, 0.5)
    else:
        zb = max(zb, 0.25)
    return zb

def zf(x):
    pente = (xb-xa/ztop-zbot)
    if x<xa:
        zf = zbot
    elif xa<=x and x<xb:
        zf = zbot + pente*(x-xa)
    elif xb<=x and x<xc:
        zf = ztop
    elif xc<=x and x<xd:
        zf = ztop - pente*(x-xc)
    else:
        zf = zbot
    return zf

def h(x):
    return max(1.2 - zb(x), 0.)

if __name__ == '__main__':

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=500,
        xa=0.,
        xb=20.,
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
        
    # build bnd condition:
    h0 = 0.95
    q0 = 0.1
    a = 0.01
    tm = 50.

    def q_bnd(t):
        if t < tm:
            return -(q0 + a*t)
        else:
            return -(q0 + 2.*a*tm -a*t)

    def h_bnd(t):
        if t < tm:
            return h0 + a*t
        else:
            return h0 + 2.*a*tm -a*t

    print("~~> create boundary condition")
    build_bnd_file(
        bnd_file="BCQ.liq", 
        ntimes=50, 
        t0=0., 
        tf=150.,
        q_funct=q_bnd, 
        verbose=False)
        
    build_bnd_file(
        bnd_file="BCH.liq", 
        ntimes=50, 
        t0=0., 
        tf=150.,
        h_funct=h_bnd, 
        verbose=False)
        
    # run init 
    print("~~> run init")
    os.system("arrakis.py bedload_weir_ini.yml")
    
    # build continuation file
    print("~~> create initial condition (continuation)")
    os.system("cp RESU/RESfin.dat .")
    os.system("mv RESfin.dat CONT.ini")
