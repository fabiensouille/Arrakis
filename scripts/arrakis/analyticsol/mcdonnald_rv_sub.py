""" 
Compute the 1D analytical solution of the Mc.Donnald test case
"""
import os
import numpy as np
import matplotlib.pylab as plt
from scipy.integrate import quad

g = 9.80665 # Gravity
Q = 20.     # Flow rate
n = 0.03    # Manning coefficient

class McDonnaldRVSubAnalyticSol():
    """
    Mc.Donnald 1d analytical solution
    Description: subcritical flow in a rectangular channel 
                 with variable width
    example 2 of : A kinetic interpretation of the section-averaged 
    Saint-Venant system for natural river hydraulics
    see also: Mac Donald 1995. Comparison of some steady state 
    Saint-Venant solvers for some test problems with analytic solutions
    """
    def __init__(self, xa=0., xb=5000., N=10000):

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
        fact1 = -64.*((x/5000.) - (1./3.))**2
        fact2 = -64.*((x/5000.) - (2./3.))**2
        L = 10. - 5.*np.exp(fact1) - 5.*np.exp(fact2)
        return L

    def depth(self, x):
        """ channel_depth """
        fact1 = -64.*((x/5000.) - (1./3.))**2
        fact2 = -64.*((x/5000.) - (2./3.))**2
        H = 1. + 0.5*np.exp(fact1) + 0.5*np.exp(fact2)
        return H

    def width_derivative(self, x):
        """ channel_width """
        fact1 = -64.*((x/5000.) - (1./3.))**2
        fact2 = -64.*((x/5000.) - (2./3.))**2
        dfact1 = -(128./5000.)*((x/5000.) - (1./3.))
        dfact2 = -(128./5000.)*((x/5000.) - (2./3.))
        dLdx = - 5.*dfact1*np.exp(fact1) - 5.*dfact2*np.exp(fact2)
        return dLdx

    def depth_derivative(self, x):
        """ channel_depth """
        fact1 = -64.*((x/5000.) - (1./3.))**2
        fact2 = -64.*((x/5000.) - (2./3.))**2
        dfact1 = -(128./5000.)*((x/5000.) - (1./3.))
        dfact2 = -(128./5000.)*((x/5000.) - (2./3.))
        dHdx = - 0.5*dfact1*np.exp(fact1) - 0.5*dfact2*np.exp(fact2)
        return dHdx

    def slope(self, x):
        """ channel_slope """
        L = self.width(x)
        H = self.depth(x)
        dL = self.width_derivative(x)
        dH = self.depth_derivative(x)
        dz = (1. - (Q**2)/(g*(L**2)*(H**3)))*dH \
           + (Q**2)*(n**2)*( ((2*H + L)**(4./3.))/((L*H)**(10./3.)) ) \
           - (Q**2)*( (dL)/(g*(L**3)*(H**2)) )
        return dz

    def bottom(self, x, xb=5000.):
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
    sol = McDonnaldRVSubAnalyticSol(xa=0., xb=5000.)
    sol()
    plt.plot(sol.x, sol.zb, label='zb')
    plt.plot(sol.x, sol.E, label='E')
    plt.legend()
    plt.show()


if __name__ == '__main__':

    # TEST ANALYTIC SOLUTION:
    test_solution()
