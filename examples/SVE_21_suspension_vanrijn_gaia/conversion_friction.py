#!/usr/bin/env python
import os
import argparse
import numpy as np
import matplotlib.pylab as plt

g = 9.81
nu = 1.e-6
kappa = 0.40

def Cf_Chezy(C):
    Cf = (2.*g)/C**2 
    return Cf

def Cf_Strickler(Ks, h):
    Cf = (2.*g)/((Ks**2)*(h**(1./3.)))
    return Cf
    
def Ks_Strickler(Cf, h):
    Ks = np.sqrt(2*g/(Cf*h**(1./3.)))
    return Ks

def Cf_ColebrookWhite(ks, h, u):
    Re = (4.*h*u)/nu
    C1 = 2.
    C2 = 2.51
    C3 = 12.
    niter_max = 50
    r = 1.
    prec = 1.e-8
    x0 = -3.6*np.log( 6.9/Re + (ks/(12.*h))**1.11 )/np.log(10.)
    niter = 1
    while (niter <= niter_max and r>=prec):
        x = -2.*C1*np.log( (C2/Re)*x0 + ks/(C3*h) )/np.log(10.)
        r = abs(x-x0)
        x0 = x
        niter += 1
    Cf = 1./x**2
    return Cf

def Cf_Haaland(ks, h, u):
    Re = (4.*h*u)/nu
    x = -3.6*np.log( 6.9/Re + (ks/(12.*h))**1.11 )/np.log(10.)
    Cf = 1./x**2
    return Cf

def Cf_Nikuradse(ks, h):
    aux = np.log(11.036*h/ks)
    aux = max(1.001, aux)
    Cf = (2.*kappa**2)/(aux**2)
    return Cf
    
def ramette(strickler=None, ks=None):
    if strickler is not None:
        return (8.3*np.sqrt(g)/(strickler))**6.
    if ks is not None:
        return 8.3*np.sqrt(g)/(ks**(1./6.))

if __name__ == '__main__':

    # Parameters gaia
    # ~~~~~~~~~~~~~~~
    h = 0.39
    d50 = 160.e-6
    alpha_ks = 1.9
    ks_hyd = 0.025
    ks_sed = alpha_ks*d50
    ks_ratio = ks_sed/ks_hyd

    print("ks_hyd   = ", ks_hyd)
    print("ks_sed   = ", ks_sed)
    print("ks_ratio = ", ks_ratio)

    # Computing ratios of Cf
    # ~~~~~~~~~~~~~~~~~~~~~~
    Cf_hyd = Cf_Nikuradse(ks_hyd, h)
    Cf_sed = Cf_Nikuradse(ks_sed, h)
    Cf_ratio = Cf_sed/Cf_hyd 

    print("Cf_hyd   = ", Cf_hyd)
    print("Cf_sed   = ", Cf_sed)
    print("Cf_ratio = ", Cf_ratio)

    # Computing ratios of Ks
    # ~~~~~~~~~~~~~~~~~~~~~~
    K_hyd = Ks_Strickler(Cf_hyd, h)
    K_sed = Ks_Strickler(Cf_sed, h)
    K_ratio = K_sed/K_hyd

    print("K_hyd    = ", K_hyd)
    print("K_sed    = ", K_sed)
    print("K_ratio  = ", K_ratio)

#    # strickler from Ramette formula:
#    K = ramette(ks=ks)
#    print("strickler from Ramette formula:", K)

#    # compute friction coef Cf from nikuradse
#    Cf_nik = Cf_Nikuradse(ks, h)
#    print("Cf_nik", Cf_nik)

#    # compute strickler coefficient from friction coef Cf
#    K_str = Ks_Strickler(Cf_nik, h)
#    print("K_str", K_str)

#    # verification
#    Cf_str = Cf_Strickler(47.16345, h)
#    print("Cf_str", Cf_str)
