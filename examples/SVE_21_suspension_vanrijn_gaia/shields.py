#!/usr/bin/env python
import os
import argparse
import numpy as np
import matplotlib.pylab as plt

Kappa = 0.4
nu = 1e-6
rho = 1000.

def settling_velocity_gaia(d50, rho_s=2650.):
    """ Compute settling velocity (Gaia) """
    if d50 <= 10**-4:
        ws = ((rho_s/rho)-1)*9.81*d50**2/(18*nu)
    elif d50 <= 10**-3:
        ws = 10*nu/d50*((1+0.01*((rho_s/rho)-1)*9.81*d50**3/nu**2)**0.5 -1)
    else:
        ws = 1.1*(((rho_s/rho)-1)*9.81*d50)**0.5
    return ws

def taus_cr(d50, rho_s=2650.):
    """ Compute critical shields """
    dadim = d50*(((rho_s-rho)/rho)*9.81/nu**2)**(1/3)
    if dadim <=4:
        ShieldsCr = 0.24*dadim**(-1)
    elif dadim <=10:
        ShieldsCr = 0.14*dadim**(-0.64)
    elif dadim <=20:
        ShieldsCr = 0.04*dadim**(-0.10)
    elif dadim <=150:
        ShieldsCr = 0.013*dadim**(0.29)
    else:
        ShieldsCr = 0.045
    return ShieldsCr

if __name__ == '__main__':

    rho_s = 2650.
    d50 = 160.e-6
    tau_star_cr = taus_cr(d50)

    print(tau_star_cr)
