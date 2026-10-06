""" 
Compute the 1D analytical solution of the traceur advection test
"""
import os
from pathlib import Path
import numpy as np
import matplotlib.pylab as plt

#------------------------
# PARAMETERS:
#------------------------
G = 9.81

#------------------------
# CLASSES:
#------------------------
class TraceurAdvection():
    """
    Thacker 1d analytical solution

    """
    def __init__(self, u0=0.1, x0=2., x1=4.5, x2=6.5, length=10., N=4001):

        # Parameters
        self.x0 = x0
        self.x1 = x1
        self.x2 = x2
        self.u0 = u0 # Advection velocity in m/s

        # 1D mesh
        self.N = N
        self.x = np.linspace(0., length, N)

        # Solution
        self.T = np.empty((N))

    def __call__(self, time):
        for i in range(self.N):

            xt0 = self.x0 + time*self.u0
            xt1 = self.x1 + time*self.u0 
            xt2 = self.x2 + time*self.u0 

            self.T[i] = np.exp(-2.5*(self.x[i]-xt0)**2)

            if (xt1 <= self.x[i] and self.x[i] <= xt2):
                self.T[i] = 1.

 
    def savetxt(self, path=''):
        filepath = os.path.join(path, 'ANALYTIC_SOL.txt') if path else 'ANALYTIC_SOL.txt'
        np.savetxt(filepath, np.c_[self.x, self.H, self.U, self.E, self.F])

    def cleantxt(self, path=''):
        file_path = Path(path) / 'ANALYTIC_SOL.txt'
        file_path.unlink(missing_ok=True)


def test_solution():
    sol = TraceurAdvection()
    sol(0.)
    plt.plot(sol.x, sol.T, label='T(0)')
    sol(10.)
    plt.plot(sol.x, sol.T, label='T(10)')
    plt.legend()
    plt.show()

#------------------------
# MAIN:
#------------------------
if __name__ == '__main__':

    # TEST ANALYTIC SOLUTION:
    test_solution()
