#!/usr/bin/env python
import os
import numpy as np
import time
import matplotlib.pylab as plt
from arrakis import *

if __name__ == '__main__':

    num_thread = [1, 2, 3, 4] 

    # config = debug
    cpu_t_0 = np.asarray([7.041422128677368, 4.142279863357544, 2.9550986289978027, 2.5217785835266113])
    init = 0.34029078483581543
    cpu_t_0 = cpu_t_0 - init
    speedup0 = cpu_t_0[0]/cpu_t_0

    # config = release
    cpu_t_1 = np.asarray([8.37524127960205, 5.153632640838623, 3.9073920249938965, 3.3347740173339844])
    init = 0.2059774398803711
    cpu_t_1 = cpu_t_1 - init
    speedup1 = cpu_t_1[0]/cpu_t_1

    # config = optim
    cpu_t_2 = np.asarray([8.7682204246521, 5.207266092300415, 3.9810235500335693, 3.5580642223358154])
    init = 0.203871488571167
    cpu_t_2 = cpu_t_2 - init
    speedup2 = cpu_t_2[0]/cpu_t_2

    # plot speed up
    set_rcparams()
    c0, c1, c2 = get_default_color_palette()
    markers = get_default_markers()
    colors = c0 + c1 + c2
    fig, ax = plt.subplots(1, 1, figsize=(4.5, 3.5))
    ax.plot(num_thread, num_thread, c='k', lw=0.5, ls='--', marker='', label="Ideal")
    ax.plot(num_thread, speedup0, c='b', lw=0.5, ls='-', marker='s', label="Arrakis")
    #ax.plot(num_thread, speedup0, c='b', lw=0.5, ls='-', marker='s', label="Arrakis (Debug)")
    #ax.plot(num_thread, speedup1, c='g', lw=0.5, ls='-', marker='s', label="Arrakis (Release)")
    #ax.plot(num_thread, speedup2, c='r', lw=0.5, ls='-', marker='s', label="Arrakis (Optim)")
    plt.xlabel("Number of cores")
    plt.ylabel("Acceleration factor")
    plt.legend()
    plt.grid()
    plt.savefig("FIG/speedup.pdf", dpi=300)
    plt.show()
