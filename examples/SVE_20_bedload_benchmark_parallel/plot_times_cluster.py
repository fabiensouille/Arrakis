#!/usr/bin/env python
import os
import numpy as np
import time
import matplotlib.pylab as plt
from arrakis import *

if __name__ == '__main__':

    num_thread = [1, 6, 12, 18, 24, 48]

    # config = release
    cpu_t_1 = np.asarray([52.42, 13.98, 9.63, 8.25, 7.57, 7.17])
    init = 0.57
    cpu_t_1 = cpu_t_1 - init
    speedup1 = cpu_t_1[0]/cpu_t_1

    # plot speed up
    set_rcparams()
    c0, c1, c2 = get_default_color_palette()
    markers = get_default_markers()
    colors = c0 + c1 + c2
    fig, ax = plt.subplots(1, 1, figsize=(4., 3.))
    ax.plot(num_thread, num_thread, c='k', lw=0.5, ls='--', marker='', label="Ideal")
    ax.plot(num_thread, speedup1, c='b', lw=0.5, ls='-', marker='s', label="Arrakis")
    plt.legend()
    plt.grid()
    plt.savefig("FIG/speedup_cluster.pdf", dpi=300)
    plt.show()
