import numpy as np
import matplotlib.pylab as plt
from arrakis import *

if __name__ == '__main__':

    # Read reference data
    Q100_Bleme = np.loadtxt("Q100_Bleme.loi", skiprows=3)
    Q100_Buech = np.loadtxt("Q100_Buech.loi", skiprows=3)
    Q100_Cote_aval = np.loadtxt("Q100_Cote_aval.loi", skiprows=3)
    
    print(Q100_Bleme)
    print(Q100_Cote_aval)

    # BC LEFT
    # *******
    # write .liq file
    f = open("BCQ1.liq", "w")
    f.write("# Liquid boundary file \n")
    f.write("T,Q \n")
    f.write("s,m^3/s \n")
    for i in range(np.shape(Q100_Buech)[0]):
        time = Q100_Buech[i, 0]
        q = Q100_Buech[i, 1]
        f.write("{:.8f} {:.8f} \n".format(time, q))
    f.close()

    # write .liq file
    f = open("BCQ2.liq", "w")
    f.write("# Liquid boundary file \n")
    f.write("T,Q \n")
    f.write("s,m^3/s \n")
    for i in range(np.shape(Q100_Bleme)[0]):
        time = Q100_Bleme[i, 0]
        q = Q100_Bleme[i, 1]
        f.write("{:.8f} {:.8f} \n".format(time, q))
    f.close()
    
    # write .liq file
    f = open("BCQ.liq", "w")
    f.write("# Liquid boundary file \n")
    f.write("T,Q \n")
    f.write("s,m^3/s \n")
    for i in range(np.shape(Q100_Bleme)[0]):
        time = Q100_Bleme[i, 0]
        q = Q100_Bleme[i, 1]+Q100_Buech[i, 1]
        f.write("{:.8f} {:.8f} \n".format(time, q))
    f.close()

    # BC RIGHT
    # ********
    # write .liq file
    f = open("BCH.liq", "w")
    f.write("# Liquid boundary file \n")
    f.write("T,H \n")
    f.write("s,m \n")
    zb_aval = 630.
    for i in range(np.shape(Q100_Cote_aval)[0]):
        time = Q100_Cote_aval[i, 0]
        h = Q100_Cote_aval[i, 1] - zb_aval
        f.write("{:.8f} {:.8f} \n".format(time, h))
    f.close()
   
    # plot
    # ****
    set_rcparams()
    c0, c1, c2 = get_default_color_palette()
    fig, ax = plt.subplots(1, 1, figsize=(5., 4.))
    # Q0
    ax.plot(Q100_Bleme[:, 0]/3600., Q100_Bleme[:, 1]+Q100_Buech[:, 1], c=c0[0], ls='-', label='$Q_0(t)$')
    #ax.plot(Q100_Bleme[:, 0]/3600., Q100_Bleme[:, 1], color=c0[3], ls='-.', label='$Q_{Bleme}$')
    #ax.plot(Q100_Buech[:, 0]/3600., Q100_Buech[:, 1], color=c0[4], ls='-.', label='$Q_{Buech}$')
    ax.set_ylabel("$Q_0(t)$ (m$^3$/s)")
    ax.set_xlabel("$t$ (h)")
    ax.set_xlim([0., 50.])
    ax.set_ylim([0., 1200.])
    ax.grid()
    ax.legend(loc=2)
    # HL
    axb = ax.twinx()
    axb.plot(Q100_Cote_aval[:, 0]/3600., Q100_Cote_aval[:, 1]-zb_aval, c=c0[1], ls='--', label='$H_L(t)$')
    axb.set_ylabel("$H_L(t)$ (m)")
    axb.set_ylim([0., 10.])
    axb.legend(loc=1)
    plt.savefig("saint_matthieu_bc.png", dpi=300)
    plt.show()
