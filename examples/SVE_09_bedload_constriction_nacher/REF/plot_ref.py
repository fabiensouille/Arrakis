import numpy as np
import matplotlib.pylab as plt
from scipy.integrate import quad
from scipy.misc import derivative

"""
Nácher-Rodríguez B, Vallés-Morán FJ, Balaguer-Beser A, Capilla MT. 
Numerical-experimental modelling of local scouring downstream of protected bridges in alluvial river beds. 
Paper presented at: E-proceedings of the 36th IAHR World Congress; 2015; The Hague, Netherlands.
"""

if __name__ == '__main__':

    # ref data
    ref_data_10 = np.loadtxt('ref10.txt')
    ref_data_30 = np.loadtxt('ref30.txt')
    ref_data_60 = np.loadtxt('ref60.txt')
    ref_data_90 = np.loadtxt('ref90.txt')

    # plot
    fig, ax = plt.subplots(1, 1, figsize=(5.,4.))
    ax.plot(ref_data_10[:, 0], ref_data_10[:, 1], label="10 min.", color="r")
    ax.plot(ref_data_30[:, 0], ref_data_30[:, 1], label="30 min.", color="b")
    ax.plot(ref_data_60[:, 0], ref_data_60[:, 1], label="60 min.", color="g")
    ax.plot(ref_data_90[:, 0], ref_data_90[:, 1], label="90 min.", color="m")
    plt.legend()
    #ax.set_ylim([zb[-1], zb[0]+h[0]])
    ax.set_xlim([0., 0.9])
    ax.set_ylabel("$z$ (m)")
    ax.set_xlabel("$x$ (m)")
    plt.show()
    plt.close()

