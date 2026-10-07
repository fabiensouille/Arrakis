#!/usr/bin/env python
import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate
from scipy.integrate import quad
from arrakis import *

G = 9.80665
Ks = 0.

if __name__ == '__main__':

    # read zb file
    data_zb = np.loadtxt('REF/cote_fond.txt')
    x0 = data_zb[:, 0]
    zb0 = data_zb[:, 1]
    print(len(x0))

    # read L file
    data_zb = np.loadtxt('REF/canal.txt')
    x1 = data_zb[:, 0]
    B1 = data_zb[:, 2] - data_zb[:, 1]
    print(len(x1))

    # interpolators
    zb_interp = interpolate.interp1d(x0, zb0)
    Lx_interp = interpolate.interp1d(x1, B1)

    # Mesh 0
    # ~~~~~~
    x = x1
    nx = len(x)
    Lx = B1
    z1 = 0.
    k = 1.
    zb = np.zeros((nx), dtype='d')
    HH = np.zeros((nx), dtype='d')
    h = np.ones((nx), dtype='d')
    q = np.ones((nx), dtype='d')
    u = np.ones((nx), dtype='d')
    for i in range(nx):
        HH[i] = z1 + ((k**2)/(2.*G))*(1./(Lx[i]**2.) - 1./(Lx[0]**2.))
        zb[i] = -HH[i]
        q[i] = h[i]*k/Lx[i]
        u[i] = k/Lx[i]

    fig, ax = plt.subplots(1, 1, figsize=(5.,4.))
    #ax.plot(x, Lx/2, color='b', ls='-', label='B')
    ax.plot(x, zb, color='b', ls='-', label='zb')
    ax.plot(x, zb+h, color='r', ls='-', label='zb+h')
    plt.legend()
    plt.show()
    plt.close()

    # write .geo file
    f = open("MESH_{}.geo".format(nx), "w")
    f.write("# Mesh file, Nx={} \n".format(nx))
    f.write("X,ZB,L \n")
    for i in range(nx):
        f.write("{:.8f} {:.8f} {:.8f} \n".format(x[i], zb[i], Lx[i]))
    f.close()
    
    # write .ini file
    f = open("INI_{}.ini".format(nx), "w")
    f.write("0., {} \n".format(nx))
    f.write(" X,H,U,0,0,0,0,0,0 \n")
    for i in range(nx):
        f.write("{:.8f} {:.8f} {:.8f} \n".format(x[i], h[i], u[i]))
    f.close()

    # Mesh 1
    # ~~~~~~
    x = x0
    nx = len(x)
    Lx = Lx_interp(x)
    z1 = 0.
    k = 1.
    zb = np.zeros((nx), dtype='d')
    HH = np.zeros((nx), dtype='d')
    h = np.ones((nx), dtype='d')
    q = np.ones((nx), dtype='d')
    u = np.ones((nx), dtype='d')
    for i in range(nx):
        HH[i] = z1 + ((k**2.)/(2.*G))*(1./(Lx[i]**2.) - 1./(Lx[0]**2.))
        zb[i] = -HH[i]
        q[i] = h[i]*k/Lx[i]
        u[i] = k/Lx[i]

    # write .geo file
    f = open("MESH_{}.geo".format(nx), "w")
    f.write("# Mesh file, Nx={} \n".format(nx))
    f.write("X,ZB,L \n")
    for i in range(nx):
        f.write("{:.8f} {:.8f} {:.8f} \n".format(x[i], zb[i], Lx[i]))
    f.close()
    
    # write .ini file
    f = open("INI_{}.ini".format(nx), "w")
    f.write("0., {} \n".format(nx))
    f.write(" X,H,U,0,0,0,0,0,0 \n")
    for i in range(nx):
        f.write("{:.8f} {:.8f} {:.8f} \n".format(x[i], h[i], u[i]))
    f.close()
