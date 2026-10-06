"""
Compute the 1D stationnary analytical solution of the bump test case
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
class Bottom():
    """
    Bottom elevation class
    """
    def __init__(self, bottom_function, x, xb=5.):
        self.function = bottom_function
        self.xb = xb

        if isinstance(x, list) or isinstance(x, np.ndarray):
            self.zb = np.empty(len(x))
            for i in range(len(x)):
                self.zb[i] = self.compute_bottom(x[i])
        else:
            self.zb = self.compute_bottom(x)

    def compute_bottom(self, x):
        if self.function == 'parabolic':
            return self.parabolic(x)
        elif self.function == 'exponential':
            return self.exponential(x)
        else:
            raise ValueError("Unknown bottom_function")

    def parabolic(self, x):
        if 4. < x < 6.:
            return 0.1 - 0.1*(x - self.xb)**2
        else:
            return 0.

    def exponential(self, x):
        #return 0.25*np.exp(-5.*(x - self.xb)**2)
        return 0.1 + 0.1*np.exp(-(x-self.xb)**2.0)

class BumpAnalyticSol():
    """
    Bumps 1d stationary analytical solution
    """
    def __init__(\
            self, flow='sub', Q=8.85, hl=1.8, length=10., width=1.,
            xb=5., bottom_function='parabolic', N=10001):

        # Physical properties
        self.flow = flow
        self.q0 = Q/width
        self.hc = (self.q0**2/G)**(1./3.)
        self.hl = hl

        # 1D mesh
        self.N = N
        self.x = np.linspace(0., length, N)

        # Bottom
        bottom = Bottom(bottom_function, self.x, xb=xb)
        self.zb = bottom.zb
        self.zm = max(self.zb)

        # Solution
        self.H = np.empty((N))
        self.E = np.empty((N))
        self.U = np.empty((N))
        self.V = np.empty((N))
        self.F = np.empty((N))

    def __call__(self):
        if self.flow == 'sub':
            self.compute_subcritical_solution()
        elif self.flow == 'cri':
            self.compute_critical_solution()
        elif self.flow == 'trans':
            self.compute_transcritical_solution()

        for i in range(self.N):
            self.E[i] = self.zb[i] + self.H[i]
            self.U[i] = self.q0 / self.H[i]
            self.V[i] = 0.0
            self.F[i] = abs(self.U[i])/np.sqrt(G*self.H[i])

    def compute_subcritical_solution(self):
        """ Subcritical solution """
        coeff = np.zeros((4), dtype='d')
        for i in range(self.N):
            polyd3_coef(coeff, self.zb[i], self.q0, self.hl, self.zb[-1])
            hroots = cardan_roots(coeff)
            hsub = max(hroots)
            # Subcritical flow for all i:
            self.H[i] = hsub

    def compute_critical_solution(self):
        """ Transcritical solution without shock """
        coeff = np.zeros((4), dtype='d')
        for i in range(self.N):
            polyd3_coef(coeff, self.zb[i], self.q0, self.hc, self.zm)
            hroots = cardan_roots(coeff)
            hsub = hroots[2]
            hsup = hroots[1]
            j = np.argmax(self.zb)
            if self.x[i] < self.x[j]:
                # Subcritical flow:
                self.H[i] = hsub
            elif self.x[i] > self.x[j]:
                # Supercritical flow:
                self.H[i] = hsup
            else:
                # Transition:
                self.H[i] = self.hc

    def compute_transcritical_solution(self):
        """ Transcritical solution with shock """
        coeff = np.zeros((4), dtype='d')
        i = self.N - 1
        supercritical = False
        hjump = False
        h1d = np.empty(self.N, dtype='d')

        while i >= 0:
            while supercritical == False:
                j = np.argmax(self.zb)
                if self.x[i] > self.x[j]:
                    polyd3_coef(coeff, self.zb[i], self.q0, self.hl, self.zb[-1])
                    hroots = np.sort(np.roots(coeff).real)
                    h1d[i] = hroots[2]
                    if h1d[i] < self.hc:
                        supercritical = True
                i = i - 1
                if i < 0: break

            if supercritical == True:
                ic = np.argmax(self.zb)
                while hjump == False:
                    polyd3_coef(coeff, self.zb[i], self.q0, self.hc, self.zm)
                    hroots = np.sort(np.roots(coeff).real)
                    h1d[ic] = hroots[1]
                    if ic > i:
                        state_L = 2*self.q0**2/h1d[ic] + G*h1d[ic]**2
                        state_R = 2*self.q0**2/h1d[ic+1] + G*h1d[ic+1]**2
                        if state_L < state_R:
                            hjump = True
                            x_shock = self.x[ic]
                            i = -1
                    ic += 1

        if not hjump:
            raise ValueError("No jump, check flow type")

        for i in range(self.N):
            if self.x[i] > x_shock:
                polyd3_coef(coeff, self.zb[i], self.q0, self.hl, self.zb[-1])
                hroots = cardan_roots(coeff)
                hsub = max(hroots)
                # Subcritical flow:
                self.H[i] = hsub

        for i in range(self.N):
            if self.x[i] <= x_shock:
                polyd3_coef(coeff, self.zb[i], self.q0, self.hc, self.zm)
                hroots = cardan_roots(coeff)
                hsub = hroots[2]
                hsup = hroots[1]
                j = np.argmax(self.zb)
                if self.x[i] < self.x[j]:
                    # Subcritical flow:
                    self.H[i] = hsub
                elif self.x[i] > self.x[j]:
                    # Supercritical flow:
                    self.H[i] = hsup
                else:
                    # Transition:
                    self.H[i] = self.hc

    def savetxt(self, path=''):
        filepath = os.path.join(path, 'ANALYTIC_SOL.txt') if path else 'ANALYTIC_SOL.txt'
        np.savetxt(filepath, np.c_[self.x, self.H, self.U, self.E, self.F])

    def cleantxt(self, path=''):
        file_path = Path(path) / 'ANALYTIC_SOL.txt'
        file_path.unlink(missing_ok=True)

#------------------------
# FUNCTIONS:
#------------------------
def polyd3_coef(coef, z, q, zs, zm):
    coef[0] = 1.0
    coef[1] = z - (q**2.0)/(2.*G*zs**2) - zs - zm
    coef[2] = 0.0
    coef[3] = (q**2.)/(2.*G)
    return 0

def cardan_roots(coef):
    # polynomial coefficients
    a1 = coef[1]
    a2 = coef[2]
    a3 = coef[3]
    # Cardan's method
    Q = (a1**2 - 3.*a2)/9.
    R = (2.*a1**3 - 9.*a1*a2 + 27.*a3)/54.
    D = Q**3 + R**2
    lamb = np.zeros((3))
    if D >= 0:
        theta = np.arccos(R/np.sqrt(Q**3))
        for k in range(3):
            lamb[k] = (-2.*np.sqrt(Q))*np.cos((theta+2.*np.pi*(k-1.))/3.)-(a1/3.)
    else:
        s = -1.*copysign(1., R)*(abs(R)+np.sqrt(R**2-Q**3))**(1./3.)
        omega = complex(-0.5, np.sqrt(3.)/2.)
        lamb[1] = s + Q/s - a1/3.
        z2 = (s*omega) + ((Q/s)*omega**2) - a1/3.
        z3 = (s*omega**2) + ((Q/s)*omega) - a1/3.
        lamb[2] = z2.real
        lamb[3] = z3.real
    return np.sort(lamb)

#------------------------
# UNIT TESTS:
#------------------------
def test_bottom_class():
    x = np.linspace(0., 20., 2001)
    bottom_functions = ['parabolic', 'exponential']
    for funct in bottom_functions:
        bottom = Bottom(funct, x, xb=10.)
        plt.plot(x, bottom.zb, label=funct)
    plt.legend()
    plt.show()

def test_solution_exp():
    bottom_function = 'exponential'
    solsub = BumpAnalyticSol(\
        flow='sub', Q=0.2, hl=0.4, bottom_function=bottom_function)
    solcri = BumpAnalyticSol(\
        flow='cri', Q=0.15, hl=0., bottom_function=bottom_function)
    soltrans = BumpAnalyticSol(\
        flow='trans', Q=0.4, hl=0.4, bottom_function=bottom_function)
    solsub()
    solcri()
    soltrans()
    plt.plot(solsub.x, solsub.zb, label='zb')
    plt.plot(solsub.x, solsub.E, label='sub')
    plt.plot(solsub.x, solcri.E, label='cri')
    plt.plot(solsub.x, soltrans.E, label='trans')
    plt.legend()
    plt.show()

def test_solution_par():
    bottom_function = 'parabolic'
    solsub = BumpAnalyticSol(\
        flow='sub', Q=0.2, hl=0.4, bottom_function=bottom_function)
    solcri = BumpAnalyticSol(\
        flow='cri', Q=0.15, hl=0., bottom_function=bottom_function)
    soltrans = BumpAnalyticSol(\
        flow='trans', Q=0.4, hl=0.4, bottom_function=bottom_function)
    solsub()
    solcri()
    soltrans()
    plt.plot(solsub.x, solsub.zb, label='zb')
    plt.plot(solsub.x, solsub.E, label='sub')
    plt.plot(solsub.x, solcri.E, label='cri')
    plt.plot(solsub.x, soltrans.E, label='trans')
    plt.legend()
    plt.show()

def debug_solution_trans():
    bottom_function = 'parabolic'
    Qin = 0.45
    HL = .38
    solsub = BumpAnalyticSol(\
        flow='sub', Q=Qin, hl=HL, bottom_function=bottom_function)
    solcri = BumpAnalyticSol(\
        flow='cri', Q=Qin, hl=HL, bottom_function=bottom_function)
    soltrans = BumpAnalyticSol(\
        flow='trans', Q=Qin, hl=HL, bottom_function=bottom_function)
    solsub()
    solcri()
    soltrans()
    plt.plot(solsub.x, solsub.zb, label='zb')
    plt.plot(solsub.x, solsub.zb+solsub.hc, label='hc')
    plt.plot(solsub.x, solsub.E, label='sub')
    plt.plot(solsub.x, solcri.E, label='cri')
    plt.plot(solsub.x, soltrans.E, label='trans')
    plt.legend()
    plt.show()

#------------------------
# MAIN:
#------------------------
if __name__ == '__main__':

    # TEST BOTTOM CLASS:
    #test_bottom_class()

    # TEST ANALYTIC SOLUTION:
    #test_solution_exp()
    test_solution_par()
    #debug_solution_trans()
