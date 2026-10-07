import os
import numpy as np
import matplotlib.pylab as plt
from arrakis import *

PLOT = True

def zb(x):
    return 5.*np.exp(-2.*(x-5.)**2.)
    
if __name__ == '__main__':

    # build mesh:
    print("~~> create mesh")
    build_geo_file(
        nx=200,
        xa=0.,
        xb=10,
        geo_file="MESH.geo",
        erodible_bed=True,
        zb_funct=zb,
        zf=-0.1,
        verbose=False)

    if PLOT:
        # plot (x,z)
        mesh_data = np.loadtxt("MESH.geo", skiprows=2)
        x = mesh_data[:,0]
        zb = mesh_data[:,1] 
        zf = mesh_data[:,2] 
        fig, ax = plt.subplots(1, 1, figsize=(5.,4.))
        ax.plot(x, zb, color='k', ls='-', lw=0.5, label='$z_b$')
        ax.plot(x, zf, color='k', ls='--', lw=0.5, label='$z_f$')
        ax.fill_between(x, zb, zf, color='peru', alpha=0.5)
        ax.fill_between(x, -0.2, zf, color='saddlebrown', alpha=0.5)
        plt.legend()
        ax.set_xlim([0., 10.])
        ax.set_ylabel("$z$ (m)")
        ax.set_xlabel("$x$ (m)")
        plt.savefig("FIG/geo_zb.png", dpi=300)
        plt.show()
        plt.close()
