import numpy as np
import matplotlib.pylab as plt
from arrakis import *

def tracers(x, k):
    " tracers initial condition, k: index for tracer"
    # tracer 1
    if k==0:
        return 0.5
    # tracer 2
    elif k==1:
        return 1.

def width(x):
    """ channel_width """
    return 1. - 0.5*np.exp(-1.*(x-5.)**2)

if __name__ == '__main__':

    # build mesh:
    # ~~~~~~~~~~~
    print("~~> create mesh")
    build_geo_file(
        nx=100,
        xa=0.,
        xb=10.,
        geo_file="MESH.geo",
        variable_width=True,
        Bx_funct=width,
        verbose=False)

    # build initial condition:
    # ~~~~~~~~~~~~~~~~~~~~~~~
    print("~~> create initial condition")
    build_ini_file(
        geo_file="MESH.geo",
        ini_file="INIC.ini",
        h=1.,
        u=0.1,
        ntrac=2,
        t_funct=tracers,
        verbose=False)
    
    # plot (x,y)
    mesh_data = np.loadtxt("MESH.geo", skiprows=2)
    ini_data = np.loadtxt("INIC.ini", skiprows=2, delimiter=',')
    x = mesh_data[:,0]
    Lx = mesh_data[:,2]
    nx = len(x)
    fig, ax = plt.subplots(1, 1, figsize=(6., 2.))
    ax.plot(x, Lx/2., color='k', ls='-', lw=0.5,  label='$\pm L/2$')
    ax.plot(x,-Lx/2., color='k', ls='-', lw=0.5, label='')
    ax.plot(x, np.zeros(nx), color='k', lw=0.5, ls=':', label='')
    ax.fill_between(x, -Lx/2., Lx/2., color='steelblue', alpha=0.25)
    plt.legend()
    ax.set_ylim([-Lx[0]/2., Lx[0]/2.])
    ax.set_xlim([0., 10.])
    ax.set_ylabel("$y$ (m)")
    ax.set_xlabel("$x$ (m)")
    plt.savefig("FIG/tracer_geo_width.png", dpi=300)
    plt.show()
    plt.close()
    

    # Build boundary condition file:
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # water level function of time
    def h(t):
        return 1.

    def u(t):
        return 0.1

    def tracers_bc(t, k):
        if k==0:
            return 0.5
        else:
            return 1.

    # write .liq file : case critical inlet
    print("~~> create boundary condition")
    build_bnd_file(
        bnd_file="BCL.liq",
        ntimes=11, 
        t0=0., 
        tf=100.,
        u_funct=u,
        ntrac=2,
        t_funct=tracers_bc,
        verbose=False)
