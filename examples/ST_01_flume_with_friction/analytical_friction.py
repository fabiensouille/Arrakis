import numpy as np
import matplotlib.pylab as plt

g = 9.80665
nu = 1.e-6
kappa = 0.41
C3 = 12. # constant in ColebrookWhite, Haaland, Nikuradse

def Cf_Chezy(C):
    Cf = (2.*g)/C**2 
    return Cf

def Cf_Strickler(Ks, Rh):
    Cf = (2.*g)/((Ks**2)*(Rh**(1./3.)))
    return Cf

def Cf_Nikuradse(ks, Rh):
    aux = max(1.001, C3*Rh/ks)
    Cf = 2./(np.log(aux)/kappa)**2.
    return Cf

def Cf_Haaland(ks, h, Rh, u):
    Re = (4.*Rh*u)/nu
    x = -3.6*np.log( 6.9/Re + (ks/(C3*Rh))**1.11 )/np.log(10.)
    Cf = 1./x**2
    return Cf

def Cf_ColebrookWhite(ks, h, Rh, u):
    Re = (4.*Rh*u)/nu
    C1 = 2.
    C2 = 2.51
    # fixed point for x:
    niter_max = 50
    r = 1.
    prec = 1.e-8
    x0 = -3.6*np.log( 6.9/Re + (ks/(C3*Rh))**1.11 )/np.log(10.)
    niter = 1
    while (niter <= niter_max and r>=prec):
        x = -2.*C1*np.log( (C2/(2.*Re))*x0 + ks/(C3*Rh) )/np.log(10.)
        r = abs(x-x0)
        x0 = x
        niter += 1
    Cf = 1./x**2
    return Cf

def Cf_Barr(ks, h, Rh, u):
    Re_star = (Rh*u)/nu
    r = ks/Rh
    aux1 = 1.1295*np.log(Re_star/1.75)/np.log(10.)
    aux2 = Re_star*(1. + (Re_star**0.52*r**0.7)/37.22)
    x = -4.*np.log( r/14.8 + aux1/aux2 )/np.log(10.)
    Cf = 1./x**2
    return Cf
    
def Cf_Bathurst(ks, h, Rh, u):
    r = ks/Rh
    x = -1.987*2.*np.log( r/5.15 )/np.log(10.)
    Cf = 1./x**2
    return Cf
    
def Cf_Machiels(ks, h, Rh, u):
    r = ks/Rh
    if r <= 0.05:
      Cf = Cf_Barr(ks, h, Rh, u)
    elif 0.05 < r <= 0.15:
      x = 1469.76*(r**3.) -382.83*(r**2.) + 9.89*r + 5.22
      Cf = 1./(2.*x)**2
    else:
      Cf = Cf_Bathurst(ks, h, Rh, u)
    return Cf

def compute_CF():
    """ Compute values of CF for different laws """
    # parameters
    S = 5./10000.
    q = 10.
    h = 5.
    u = q/h
    B = 10.
    Q = q*B
    Rh = B*h/(2*h + B)
    print('U = ', u)
    print('Q = ', Q)
    print('h = ', h)
    print('Rh = ', Rh)
    print('Fr = ', u/np.sqrt(g*h))

    # analytical friction coefficient
    Cf = 2.*g*Rh*S/(u**2)
    print("Cf                  = ", Cf)

    # friction laws
    Cf_law = Cf_Chezy(56.57)
    print("===================== ")
    print("Cf (Chezy)          = ", Cf_law)
    print("Residual            = ", abs(Cf-Cf_law))
    
    Cf_law = Cf_Strickler(48.56, Rh)
    print("===================== ")
    print("Cf (Strickler)      = ", Cf_law)
    print("Residual            = ", abs(Cf-Cf_law))

    Cf_law = Cf_Nikuradse(0.01822, Rh)
    print("===================== ")
    print("Cf (Nikuradse)      = ", Cf_law)
    print("Residual            = ", abs(Cf-Cf_law))

    Cf_law = Cf_Haaland(0.01905, h, Rh, u)
    print("===================== ")
    print("Cf (Haaland)        = ", Cf_law)
    print("Residual            = ", abs(Cf-Cf_law))

    Cf_law = Cf_ColebrookWhite(0.0192, h, Rh, u)
    print("===================== ")
    print("Cf (ColebrookWhite) = ", Cf_law)
    print("Residual            = ", abs(Cf-Cf_law))

def plot_cf():
    # parameters
    S = 5./10000.
    q = 1.
    h = 1.
    u = q/h
    B = 1.
    Rh = B*h/(2*h + B)
    
    npts = 10000
    ks = np.linspace(0.001, 1., npts)
    
    # compute Cf
    Cf_Nik = np.empty(npts)
    Cf_Col = np.empty(npts)
    Cf_Haa = np.empty(npts)
    Cf_Bar = np.empty(npts)
    Cf_Bat = np.empty(npts)
    Cf_Mac = np.empty(npts)
    
    for i in range(npts):
        k = ks[i]
        Cf_Nik[i] = Cf_Nikuradse(k, Rh)
        Cf_Haa[i] = Cf_Haaland(k, h, Rh, u)
        Cf_Col[i] = Cf_ColebrookWhite(k, h, Rh, u)
        Cf_Bar[i] = Cf_Barr(k, h, Rh, u)
        Cf_Bat[i] = Cf_Bathurst(k, h, Rh, u)
        Cf_Mac[i] = Cf_Machiels(k, h, Rh, u)
    
    # plot
    fig, ax = plt.subplots(1, 1, figsize=(6.,5.))
    ax.plot(ks[:]/h, 1./(2.*np.sqrt(Cf_Nik[:])), color='r', lw=1.5, ls='-', label='Nikuradse')
    ax.plot(ks[:]/h, 1./(2.*np.sqrt(Cf_Haa[:])), color='b', lw=1.5, ls='-', label='Haaland')
    ax.plot(ks[:]/h, 1./(2.*np.sqrt(Cf_Col[:])), color='k', lw=1.5, ls='-', label='Colebrook')
    ax.plot(ks[:]/h, 1./(2.*np.sqrt(Cf_Bar[:])), color='g', lw=0.5, ls='--', label='Barr')
    ax.plot(ks[:]/h, 1./(2.*np.sqrt(Cf_Bat[:])), color='g', lw=0.5, ls='-.', label='Bathurst')
    ax.plot(ks[:]/h, 1./(2.*np.sqrt(Cf_Mac[:])), color='g', lw=1.5, ls='-', label='Machiels')
    plt.legend()
    ax.set_xscale("log")
    #ax.set_ylim([0. , 1.])
    #ax.set_xlim([0., 5000.])
    ax.set_ylabel("$1/2 \sqrt{C_f}$")
    ax.set_xlabel("$k_s/h$")
    plt.grid()
    plt.savefig("FIG/analytical_friction.png", dpi=300)
    plt.show()
    plt.close()

if __name__ == "__main__":

    # Test CF values for vnv 
    # ~~~~~~~~~~~~~~~~~~~~~~
    compute_CF()
    
    # Plot
    # ~~~~
    plot_cf()
    

