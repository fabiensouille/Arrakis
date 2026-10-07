#!/usr/bin/env python
import os
import numpy as np
import time
import matplotlib.pylab as plt
from arrakis import *

if __name__ == '__main__':

    num_thread = [1, 2, 3, 4]

    # config = debug
    cpu_t_0 = np.asarray([78.25232005119324, 44.13425087928772, 38.09726119041443, 32.401355504989624])
    init = 1.31
    cpu_t_0 = cpu_t_0 - init
    speedup0 = cpu_t_0[0]/cpu_t_0

    # config = release
    cpu_t_1 = np.asarray([28.61041831970215, 15.395263433456421, 12.634739398956299, 10.410573244094849])
    init = 0.21884894371032715
    cpu_t_1 = cpu_t_1 - init
    speedup1 = cpu_t_1[0]/cpu_t_1

    # config = optim
    cpu_t_2 = np.asarray([25.870559215545654, 14.649511337280273, 11.78730034828186, 9.855145931243896])
    init = 0.21889615058898926
    cpu_t_2 = cpu_t_2 - init
    speedup2 = cpu_t_2[0]/cpu_t_2

    # plot speed up
    set_rcparams()
    c0, c1, c2 = get_default_color_palette()
    markers = get_default_markers()
    colors = c0 + c1 + c2
    fig, ax = plt.subplots(1, 1, figsize=(4., 3.))
    ax.plot(num_thread, num_thread, c='k', lw=0.5, ls='--', marker='', label="Ideal")
    ax.plot(num_thread, speedup1, c='b', lw=0.5, ls='-', marker='s', label="Arrakis")
    #ax.plot(num_thread, speedup0, c='b', lw=0.5, ls='-', marker='s', label="Arrakis (Debug)")
    #ax.plot(num_thread, speedup1, c='g', lw=0.5, ls='-', marker='s', label="Arrakis (Release)")
    #ax.plot(num_thread, speedup2, c='r', lw=0.5, ls='-', marker='s', label="Arrakis (Optim)")
    plt.legend()
    plt.grid()
    plt.savefig("FIG/speedup.pdf", dpi=300)
    plt.show()
