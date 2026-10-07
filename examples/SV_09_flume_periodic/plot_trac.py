import os
import matplotlib.pylab as plt
from arrakis import *

if __name__ == '__main__':

    # get res
    i = -1
    time, res, _ = load_res(resname="RES", record=i)
    time, restr1, _ = load_res(resname="RES_TRAC0001_", record=i)
    time, restr2, _ = load_res(resname="RES_TRAC0002_", record=i)

    # plot
    set_rcparams()
    fig, ax = plt.subplots(1, 1, figsize=(12, 4))
    x = restr1[:,0]
    T1 = restr1[:,1]
    T2 = restr2[:,1]
    ax.plot(x, T1, label='T1', color='b', marker='+')
    ax.plot(x, T2, label='T2', color='r', marker='+')
    plt.legend()
    plt.grid()
    plt.show()

