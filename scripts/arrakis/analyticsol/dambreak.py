"""
Compute the 1D analytical solution of the dambreak test case
"""
import os
from pathlib import Path
import numpy as np
import matplotlib.pylab as plt

#------------------------
# PARAMETERS:
#------------------------
g = 9.80665

#------------------------
# CLASSES:
#------------------------
class DambreakAnalyticSol():
    """
    Dambreak 1d analytical solution
    """
    def __init__(self, h_l=1., h_r=0.2, x_d=1, length=10., N=8001):

        # Properties
        self.x_d = x_d
        self.h_l = h_l
        self.h_r = h_r

        # Analytic type (Ritter or Stoker)
        if self.h_r == 0.:
            self.analytic_type = 1
        else:
            self.analytic_type = 2

        # 1D mesh
        self.N = N
        self.x = np.linspace(0., length, N)

        # Bottom
        self.zb = np.zeros((N))

        # Solution
        self.H = np.empty((N))
        self.E = np.empty((N))
        self.U = np.empty((N))

    def __call__(self, time):

        #==================
        #CASE 1: dry bottom
        #==================
        if self.analytic_type == 1:
            x_a = self.x_d -   time*np.sqrt(g*self.h_l)
            x_b = self.x_d + 2*time*np.sqrt(g*self.h_l)
            for i in range(self.N):
                if self.x[i] <= x_a:
                    self.H[i] = h_l
                    self.U[i] = 0.
                elif x_a <= self.x[i] and self.x[i] <= x_b:
                    self.H[i] = (4./(9.*g))*(np.sqrt(g*self.h_l)-((self.x[i]-self.x_d)/(2.*time)))**2
                    self.U[i] = (2./3.)*((self.x[i]-self.x_d)/time + np.sqrt(g*self.h_l))
                elif x_b <= self.x[i]:
                    self.H[i] = 0.
                    self.U[i] = 0.

        #==================
        #CASE 2: wet bottom
        #==================
        elif self.analytic_type == 2:

            # cm calculus:
            coeff  = np.zeros( (7), dtype='d')
            rrange = np.zeros( (2), dtype='d')
            rrange[0] = np.sqrt(g*self.h_l)
            rrange[1] = np.sqrt(g*self.h_r)
            coeff[0] = 1.0
            coeff[1] = 0.0
            coeff[2] = -9.*g*self.h_r
            coeff[3] = 16.*g*self.h_r*np.sqrt(g*self.h_l)
            coeff[4] = g**2*self.h_r*(-8.*self.h_l-self.h_r)
            coeff[5] = 0.0
            coeff[6] = (g*self.h_r)**3

            cm_r = np.roots(coeff)
            cm = max(cm_r.real)
            for i in range(len(cm_r)):
                if cm_r[i].imag == 0. and rrange[1]<= cm_r[i].real <=rrange[0]:
                    cm = cm_r[i].real
            self.check_root(rrange, cm)

            # domains 
            x_a = self.x_d - time*np.sqrt(g*self.h_l)
            x_b = self.x_d + time*(2.*np.sqrt(g*self.h_l)-3.*cm)
            x_c = self.x_d + time*((2.*cm**2*(np.sqrt(g*self.h_l)-cm))/(cm**2-g*self.h_r))

            # Analytic solution
            for i in range(self.N):
                if self.x[i] <= x_a:
                    self.H[i] = self.h_l
                    self.U[i] = 0.
                elif x_a <= self.x[i] and self.x[i] <= x_b:
                    self.H[i] = (4./(9.*g))*(np.sqrt(g*self.h_l)-((self.x[i]-self.x_d)/(2.*time)))**2
                    self.U[i] = (2./3.)*((self.x[i]-self.x_d)/time + np.sqrt(g*self.h_l))
                elif x_b <= self.x[i] and self.x[i] <= x_c:
                    self.H[i] = cm**2/g
                    self.U[i] = 2.*(np.sqrt(g*self.h_l)-cm)
                elif x_c <= self.x[i]:
                    self.H[i] = self.h_r
                    self.U[i] = 0.
            
        for i in range(self.N):
            self.E[i] = self.H[i]

    def check_root(self, rrange, cm):
        if not rrange[1]<= cm <=rrange[0]:
            raise ValueError("Failed to compute analytic solution, check roots")

    def savetxt(self, path=''):
        filepath = os.path.join(path, 'ANALYTIC_SOL.txt') if path else 'ANALYTIC_SOL.txt'
        np.savetxt(filepath, np.c_[self.x, self.H, self.U, self.E])

    def cleantxt(self, path=''):
        file_path = Path(path) / 'ANALYTIC_SOL.txt'
        file_path.unlink(missing_ok=True)


def test_solution():
    sol = DambreakAnalyticSol(h_l=1., h_r=0.2, x_d=5, length=10.)
    sol(0.5)
    plt.plot(sol.x, sol.E, label='E(0.5)')
    plt.legend()
    plt.show()

#------------------------
# MAIN:
#------------------------
if __name__ == '__main__':

    # TEST ANALYTIC SOLUTION:
    test_solution()
