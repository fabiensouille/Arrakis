import numpy as np
import matplotlib.pylab as plt
from scipy.integrate import quad
from scipy.misc import derivative

"""
Nácher-Rodríguez B, Vallés-Morán FJ, Balaguer-Beser A, Capilla MT. 
Numerical-experimental modelling of local scouring downstream of protected bridges in alluvial river beds. 
Paper presented at: E-proceedings of the 36th IAHR World Congress; 2015; The Hague, Netherlands.
"""

# Global parameters
g = 9.80665 # Gravity

def set_1d_fv_mesh(nx=100, xa=0., xb=10., density=None):
    """
    Define 1D finite volume mesh
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    The mesh has nx finite volume cells along x-axis
    Returns : x[1:nx];
    Remark: x[0]/x[nx+1] are the left/right ghost cells

    Convention:
    0    A    1        i-1        i        i+1        n    B   n+1
    o----|----o-- ... --o----|----o----|----o-- ... --o----|----o
     x0  I0   x1       xi-1  Ii-1  xi  Ii   xi+1      xn   In   xn+1
     <-------> <------->           <--------> dx[i]    <-------->
        dx[0]    dx[1]         <-------->                 dx[n]
                                 ci[i] : cell size
    Ghost cells:
    o--|--o--    ...      --o--|--o
    0  A  1                nx  B nx+1
    """
    mean_dx = (xb-xa)/nx
    x = np.zeros(nx+2)
    dx = np.zeros(nx)

    # Constant cell size
    if density is None:
        for i in range(0, nx+2):
            x[i] = 0.5*mean_dx + (i-1)*mean_dx

    return x[1:nx+1]

def width(x):
    """ channel_width """
    width_entry = 0.064
    width_middle = 0.038
    pente = (width_entry-width_middle)/0.052
    if x < 0.3705:
        L = width_entry
    elif x >= 0.3705 and x <= 0.4225:
        L = width_entry - pente*(x-0.3705)
    elif x >= 0.4225 and x <= 0.4775:
        L = width_middle
    elif x >= 0.4775 and x <= 0.5295:
        L = width_middle + pente*(x-0.4775)
    else:
        L = width_entry
    return L

def bottom(x):
    zb = 0.072
    zf = 0.04
    return zb, zf

if __name__ == '__main__':

    # Uniform mesh
    #~~~~~~~~~~~~~
    # nodes 
    x = set_1d_fv_mesh(nx=200, xa=0., xb=0.9, density=None)
    nx = len(x)
    print(x)

    # width
    Lx = np.zeros(nx)
    for i in range(nx):
        Lx[i] = width(x[i])
    print("Lx = ", Lx)
    fig, ax = plt.subplots(1, 1, figsize=(6.,3.))
    ax.plot(x, Lx/2., color='k', ls='-', lw=0.5,  label='$\pm L/2$')
    ax.plot(x,-Lx/2., color='k', ls='-', lw=0.5, label='')
    ax.plot(x, np.zeros(nx), color='k', lw=0.5, ls=':', label='')
    ax.fill_between(x, -Lx/2., Lx/2., color='steelblue', alpha=0.25)
    plt.legend()
    ax.set_ylim([-Lx[0]/2., Lx[0]/2.])
    ax.set_xlim([0., 0.9])
    ax.set_ylabel("$y$ (m)")
    ax.set_xlabel("$x$ (m)")
    plt.savefig("FIG/geo_width_{}.png".format(nx), dpi=300)
    plt.show()
    plt.close()

    # bottom
    zb = np.zeros(nx)
    zf = np.zeros(nx)
    for i in range(nx):
        zb[i], zf[i] = bottom(x[i])
    
    # depth / velocity
    h = 0.066*np.ones(nx)
    Q = 8.8333e-4
    q = np.zeros(nx)
    u = np.zeros(nx)
    for i in range(nx):
        q[i] = Q/Lx[i]
        u[i] = q[i]/h[i]
    
    print(q[0])
    print(u[0])
    
    fig, ax = plt.subplots(1, 1, figsize=(5.,4.))
    ax.plot(x, zb+h, color='b', ls='-', lw=0.5, label='$h_0$')
    ax.plot(x, zb, color='k', ls='-', lw=0.5, label='$z_b$')
    ax.fill_between(x, np.min(zb), zb, color='grey', alpha=0.5)
    ax.fill_between(x, zb, zb+h, color='steelblue', alpha=0.25)
    plt.legend()
    ax.set_ylim([zb[-1], zb[0]+h[0]])
    ax.set_xlim([0., 0.9])
    ax.set_ylabel("$z$ (m)")
    ax.set_xlabel("$x$ (m)")
    plt.savefig("FIG/geo_zb_{}.png".format(nx), dpi=300)
    plt.show()
    plt.close()

    # write .geo file
    f = open("MESH_{}.geo".format(nx), "w")
    f.write("# Mesh file, Nx={} \n".format(nx))
    f.write("X,ZB,L \n")
    for i in range(nx):
        f.write("{:.8f} {:.8f} {:.8f} \n".format(x[i], zb[i], Lx[i]))
    f.close()
    
    f = open("MESH_BEDLOAD_{}.geo".format(nx), "w")
    f.write("# Mesh file, Nx={} \n".format(nx))
    f.write("X,ZB,ZF,L \n")
    for i in range(nx):
        f.write("{:.8f} {:.8f} {:.8f} {:.8f} \n".format(x[i], zb[i], zf[i], Lx[i]))
    f.close()
    
    # write .ini file
    f = open("INI_{}.ini".format(nx), "w")
    f.write("0., {} \n".format(nx))
    f.write(" X,H,U,0,0,0,0,0,0 \n")
    for i in range(nx):
        f.write("{:.8f} {:.8f} {:.8f} \n".format(x[i], h[i], u[i]))
    f.close()
