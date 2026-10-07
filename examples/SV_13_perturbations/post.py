import os
import numpy as np
import matplotlib.pylab as plt

PLOT  = True
CLEANING = False

if __name__ == '__main__':

    # post processing:
    if PLOT:
        # get res
        file_ini = 'RESU/RES0000.dat'
        file_fin = 'RESU/RESfin.dat'
        file_ref = 'REF.dat'
        ini = np.loadtxt(file_ini, delimiter=',', skiprows=1)
        res = np.loadtxt(file_fin, delimiter=',', skiprows=1)
        ref = np.loadtxt(file_ref, delimiter=',', skiprows=1)

        # extract values
        h0 = ini[:, 2]
        u0 = ini[:, 3]
        xr = ref[:, 0]
        zr= ref[:, 1]
        hr = ref[:, 2]
        ur = ref[:, 3]
        x = res[:, 0]
        zb= res[:, 1]
        h = res[:, 2]
        u = res[:, 3]
        fr0= abs(u0)/np.sqrt(9.81*h0)
        fr= abs(u)/np.sqrt(9.81*h)
          
        # plot
        fig, axes = plt.subplots(2, 1, figsize=(6, 6))
        ax1 = axes[0]
        ax2 = axes[1]
        ax1.plot(xr, zr+hr, label='h (Reference)', c='k', ls=':', lw=1.)
        ax1.plot(x, zb+h, label='h', lw=1., marker='o', markersize=1)
        ax2.plot(xr, hr*ur, label='q (Reference)', c='k', ls=':', lw=1.)
        ax2.plot(x, h*u, label='q', lw=1., marker='s', markersize=1)

        plt.legend()
        plt.savefig("FIG/plot_res.png", dpi=300)
        plt.show()

    # cleaning:
    if CLEANING:
        os.system("rm RESU/*")
