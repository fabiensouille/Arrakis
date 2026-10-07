#!/usr/bin/env python
import numpy as np
import matplotlib.pylab as plt
from scipy import interpolate
from scipy.integrate import quad
from arrakis import *

if __name__ == '__main__':

    # read zb file
    data_zb = np.loadtxt('REF/cote_fond.txt')
    x0 = data_zb[:, 0]
    zb0 = data_zb[:, 1]
    print(len(x0))

    fig, ax = plt.subplots(1, 1, figsize=(5.,4.))
    ax.plot(x0, zb0, color='k', ls='-', label='zb')
    plt.legend()
    plt.show()
    plt.close()
    
    # read L file
    data_zb = np.loadtxt('REF/canal.txt')
    x1 = data_zb[:, 0]
    B1 = data_zb[:, 2] - data_zb[:, 1]
    print(len(x1))

    fig, ax = plt.subplots(1, 1, figsize=(5.,4.))
    ax.plot(x1, B1/2, color='b', ls='-', label='+B/2')
    ax.plot(x1,-B1/2, color='r', ls='-', label='-B/2')
    plt.legend()
    plt.show()
    plt.close()

    # interpolators
    zb_interp = interpolate.interp1d(x0, zb0)
    Lx_interp = interpolate.interp1d(x1, B1)

    # Mesh 0
    # ~~~~~~
    x = x1
    nx = len(x)
    zb = zb_interp(x)
    Lx = B1
    h = np.ones((nx), dtype='d')
    q = np.zeros((nx), dtype='d')
    u = np.zeros((nx), dtype='d')
    for i in range(nx):
        h[i] = 12. - zb[i]

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
    zb = zb0
    Lx = Lx_interp(x)
    h = np.ones((nx), dtype='d')
    q = np.zeros((nx), dtype='d')
    u = np.zeros((nx), dtype='d')
    for i in range(nx):
        h[i] = 12. - zb[i]

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
