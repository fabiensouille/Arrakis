""" 
Compute the 1D analytical solution of the Mc.Donnald test case
"""
import os
import numpy as np
import matplotlib.pylab as plt
from scipy.integrate import quad
try:
    from scipy.misc import derivative
except ImportError:
    import numdifftools as nd
    def derivative(func, x0, dx=1.0):
        return nd.Derivative(func, step=dx)(x0)

g = 9.80665 # Gravity
Q = 20.     # Flow rate
n = 0.03    # Manning coefficient

class McDonnaldRVTransAnalyticSol():
    """
    Mc.Donnald 1d analytical solution
    Description: subcritical flow in a rectangular channel 
                 with variable width
    example 2 of : A kinetic interpretation of the section-averaged 
    Saint-Venant system for natural river hydraulics
    see also: Mac Donald 1995. Comparison of some steady state 
    Saint-Venant solvers for some test problems with analytic solutions
    """
    def __init__(self, xa=0., xb=1000., N=10000):

        # 1D mesh
        self.N = N
        self.x = np.linspace(xa, xb, N)

        # Bottom
        self.zb = np.zeros((N))

        # width
        self.Lx = np.zeros((N))
        for i in range(N):
            self.Lx[i] = self.width(self.x[i])

        # Solution
        self.QL = Q*np.ones((N)) # m3/s (QL=Q*L)
        self.H = np.empty((N))
        self.E = np.empty((N))
        self.U = np.empty((N))
        self.Q = np.empty((N))   # m2/s

    def width(self, x):
        """ channel_width """
        xa = x/1000.
        L = 10. - 64.*( xa**2 -2.*xa**3 + xa**4 )
        return L

    def depth(self, x):
        """ channel_depth """
        xa = x/1000.
        if x<500.:
           H = -(1./40.) + (1./(1. + 2.*(xa-0.5)**2))
        else:
            a0 = 1.5
            ak = [-0.230680, 0.248267, -0.228271]
            H = a0*np.exp((x/4000.)-(1./4.))
            for k in range(3):
              H += ak[k]*np.exp(15.*(k+1) - 30.*(k+1)*xa)
        return H

    def slope(self, x):
        """ channel_slope """
        L = self.width(x)
        H = self.depth(x)
        Q2 = Q**2
        L2 = L**2
        L3 = L**3
        H2 = H**2
        H3 = H**3
        dL = derivative(self.width, x, dx=0.00001)
        dH = derivative(self.depth, x, dx=0.00001)
        u2 = Q2/(H2*L2)
        aux = ((2./L) + (1./H))
        Sf = (u2*n**2)*( np.sign(aux)*(np.abs(aux))**(4/3) )
        S0 = Sf + (1. - u2/(g*H))*dH - (u2/g)*(dL/L)
        return S0

    def bottom(self, x, xb=1000.):
        """ channel_bottom """
        zb = quad(func=self.slope, a=x, b=xb)
        zb = zb[0]
        return zb

    def __call__(self):

        for i in range(self.N):
            self.zb[i] = self.bottom(self.x[i])
            self.H[i] = self.depth(self.x[i])
            self.U[i] = self.QL[i]/(self.Lx[i]*self.H[i])
            self.Q[i] = self.QL[i]/(self.Lx[i])
            self.E[i] = self.zb[i] + self.H[i]

    def savetxt(self, filename='ANALYTIC_SOL.txt'):
        np.savetxt(filename,\
            np.c_[self.x, self.zb, self.H, self.U, self.Lx, self.Q, self.QL])

def test_solution():
    sol = McDonnaldRVTransAnalyticSol(xa=0., xb=1000.)
    sol()
    plt.plot(sol.x, sol.zb, label='zb')
    plt.plot(sol.x, sol.E, label='E')
    plt.legend()
    plt.show()


if __name__ == '__main__':

    # TEST ANALYTIC SOLUTION:
    test_solution()
