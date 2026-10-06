""" 
Compute the 1D analytical solution of the thacker test case
"""
import os
from pathlib import Path
import numpy as np
import matplotlib.pylab as plt

#------------------------
# PARAMETERS:
#------------------------
G = 9.80665

#------------------------
# CLASSES:
#------------------------
class ThackerAnalyticSol():
    """
    Thacker 1d analytical solution

    """
    def __init__(self, a=1, h0=0.5, length=10., N=8001):

        # Properties
        self.a = a
        self.h0 = h0
        self.length = length
        self.omega = np.sqrt(2.*G*h0)/a
        self.period = 2.*np.pi/self.omega

        # 1D mesh
        self.N = N
        self.x = np.linspace(0., length, N)

        # Bottom
        self.zb = np.zeros((N))

        # Solution
        self.H = np.empty((N))
        self.E = np.empty((N))
        self.U = np.empty((N))
        self.F = np.zeros((N))

    def __call__(self, time):
        for i in range(self.N):
            self.zb[i] = self.h0*((self.x[i]-0.5*self.length)**2/self.a**2 - 1.)
            
            x1 = -0.5*np.cos(self.omega*time) - self.a + self.length/2.
            x2 = -0.5*np.cos(self.omega*time) + self.a + self.length/2.

            if (x1 <= self.x[i] and self.x[i] <= x2):
                temp1 = (self.x[i] - self.length/2.)/self.a
                temp2 = np.cos(self.omega*time)/(2.*self.a)
                self.H[i] = -self.h0*( (temp1 + temp2)**2 - 1. )
                self.U[i] = 0.5*self.omega*np.sin(self.omega*time)
            else:
                self.H[i] = 0.
                self.U[i] = 0.

            self.E[i] = self.H[i] + self.zb[i]
 
    def savetxt(self, path=''):
        filepath = os.path.join(path, 'ANALYTIC_SOL.txt') if path else 'ANALYTIC_SOL.txt'
        np.savetxt(filepath, np.c_[self.x, self.H, self.U, self.E, self.F])

    def cleantxt(self, path=''):
        file_path = Path(path) / 'ANALYTIC_SOL.txt'
        file_path.unlink(missing_ok=True)


def test_solution():
    sol = ThackerAnalyticSol(a=3, h0=.5, length=10.)
    sol(0.)
    print(sol.period)
    print(sol.period/2.)
    plt.plot(sol.x, sol.zb, label='zb')
    plt.plot(sol.x, sol.E, label='E(0)')
    sol(sol.period/2.)
    plt.plot(sol.x, sol.E, label='E(T/2)')
    plt.legend()
    plt.show()

#------------------------
# MAIN:
#------------------------
if __name__ == '__main__':

    # TEST ANALYTIC SOLUTION:
    test_solution()
