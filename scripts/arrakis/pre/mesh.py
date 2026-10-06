import numpy as np

class Mesh1D():
    """
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Mesh 1D
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    @param nx (int): number of cells (default: 100)
    @param xa (float): absc. of left boundary interface (default: 0.)
    @param xb (float): absc. of right boundary interface (default: 10.)
    @param density[0:nx+1] (float): array of density (default: None)
                                    used to define cell size
    @attribute nx (float): number of cells
    @attribute x [0:nx+2] (float): abscissa of cells (center)
    @attribute dx[0:nx+1] (float): space between cell centers
    @attribute xa (float): absc. of left boundary interface
    @attribute xb (float): absc. of right boundary interface
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Convention:

    0   xA    1        i-1        i        i+1        n   xB   n+1
    o----|----o-- ... --o----|----o----|----o-- ... --o----|----o
     x0  I0   x1       xi-1  Ii-1  xi  Ii   xi+1      xn   In   xn+1
     <-------> <------->           <--------> dx[i]   <--------->
        dx[0]    dx[1]        <-------->                 dx[n]
                                 ci[i] : cell size

    Ghost cells: x[0] and x[nx+1]
    o--|--o--    ...      --o--|--o
    0  A  1                nx  B nx+1

    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """
    def __init__(self, nx=100, xa=0., xb=10., density=None, verbose=False):
        """
        Initialize mesh
        """
        const_dx = (xb-xa)/nx
        self.nx = nx
        self.xa = xa
        self.xb = xb
        self.x = np.zeros(nx+2)
        self.dx = np.zeros(nx+1)

        # Constant cell size
        # ~~~~~~~~~~~~~~~~~~
        for i in range(0, self.nx+2):
            self.x[i] = self.xa + 0.5*const_dx + (i-1)*const_dx
            if i>0:
                self.dx[i-1] = self.x[i] - self.x[i-1]

        # Variable cell size (keep same nx)
        # ~~~~~~~~~~~~~~~~~~
        if density is not None:
            # deform with density function
            eps_density = 1.e-3 # to avoid huge dx
            for i in range(0, self.nx+1):
                self.dx[i] /= max(eps_density, density(self.x[i]))

            # compute new xb and scaling factor
            self.x[0] = self.xa - 0.5*self.dx[0]
            for i in range(0, nx+1):
                self.x[i+1] = self.x[i] + self.dx[i]
            newxb = self.x[-1] - 0.5*self.dx[-1]
            scale = (self.xb-self.xa)/(newxb-self.xa)

            # scale
            for i in range(0, self.nx+2):
                self.x[i] = self.xa + (self.x[i]-self.xa)*scale

            # final checks
            assert self.x[1]  > self.xa
            assert self.x[nx] < self.xb

        # verbose
        if verbose:
            print("x = :", self.x[1:nx+1])

def refine_mesh(x, dx, option):
    """
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Refine mesh
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    @param x[1:nx] (float): abscissa of cell
    @param dx (float): target cell size
    @param option (int): 1: keep old mesh nodes
                         2: uniform mesh (don't keep old mesh nodes)
    ! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    """
    nx = len(x)

    # opt 1: keep old mesh nodes
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~
    if option==1:
        # define new dx by section
        nx_sec = np.empty(nx-1, dtype=int)
        dx_sec = np.empty(nx-1, dtype='d')

        for i in range(0, nx-1):
            nx_sec[i] = max(1, int((x[i+1]-x[i])/dx))
            dx_sec[i] = (x[i+1]-x[i])/nx_sec[i]

        # new cells
        nx_tot = np.sum(nx_sec)
        x_new = np.empty(nx_tot)

        for i in range(0, nx-1):
            for j in range(0, nx_sec[i]):
                if i==0 and j==0:
                    x_new[0] = x[0]
                else:
                    idx = np.sum(nx_sec[0:i]) + j
                    x_new[idx] = x[i] + j*dx_sec[i]

    # opt2: uniform mesh (don't keep old mesh nodes)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    elif option==2:
        np_new = int((x[-1]-x[0])/dx)
        dx_new = (x[-1]-x[0])/np_new
        x_new = np.arange(x[0], x[-1], dx_new)

    else:
        raise ValueError("option must be 1 or 2")
    
    return x_new