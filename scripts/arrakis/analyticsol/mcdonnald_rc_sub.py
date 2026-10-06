"""
Compute the 1D analytical solution of the Mc.Donnald test case
"""
import os
import numpy as np
import matplotlib.pylab as plt
from scipy.integrate import quad

g = 9.80665

class McDonnaldRCSubAnalyticSol():
    """
    Mc.Donnald 1d analytical solution
    Description: subcritical flow in a rectangular channel with constant width
    example 1 of : A kinetic interpretation of the section-averaged 
    Saint-Venant system for natural river hydraulics
    """
    def __init__(self, xa=0., xb=1000., N=10001):

        # 1D mesh
        self.N = N
        self.x = np.linspace(xa, xb, N)

        # Bottom
        self.zb = np.zeros((N))

        # width
        self.Lx = 10.*np.ones(N)

        # Solution
        self.QL = 20.*np.ones(N)
        self.H = np.empty((N))
        self.E = np.empty((N))
        self.U = np.empty((N))
        self.Q = np.empty((N))

    def depth(self, x):
        """ channel_depth """
        fact = -16.*((x/1000.) - (1./2.))**2
        H = ((4./g)**(1./3.))*( 1. + 0.5*np.exp(fact) )
        return H

    def depth_derivative(self, x):
        """ channel_depth """
        fact = -16.*((x/1000.) - (1./2.))**2
        dfact = -(32./1000.)*((x/1000.) - (1./2.))
        dHdx = ((4./g)**(1./3.))*0.5*dfact*np.exp(fact)
        return dHdx

    def slope(self, x):
        """ channel_slope """
        H = self.depth(x)
        dH = self.depth_derivative(x)
        dz = (1. - (4./(g*H**3)))*dH \
           + 0.36*( ((10. + 2*H)**(4./3.))/((10.*H)**(10./3.)) )
        return dz

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
    sol = McDonnaldRCSubAnalyticSol(xa=0., xb=1000.)
    sol()
    plt.plot(sol.x, sol.zb, label='zb')
    plt.plot(sol.x, sol.E, label='E')
    plt.legend()
    plt.show()


if __name__ == '__main__':

    # TEST ANALYTIC SOLUTION:
    test_solution()
