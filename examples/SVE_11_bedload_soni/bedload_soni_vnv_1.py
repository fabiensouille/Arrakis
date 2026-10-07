#!/usr/bin/env python
import os
import argparse
import numpy as np
import matplotlib.pylab as plt
from arrakis import *

class bedload_soni_1(Study):

    def pre(self):
        self.set_dir(file=__file__)
        self.name = "bedload_soni_1"
        self.cases = []

        # ACU-SWE / Nickalls' roots / Grass formula
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "bedload_soni.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'ACU-SVE'
        case_inputs['NUMERICAL PARAMETERS']['ADVECTION SCHEME'] = 3
        self.cases.append(case_inputs)

        # ROE-SWE / Cardano's roots / Grass formula
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "bedload_soni.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'ROE-SVE'
        case_inputs['NUMERICAL PARAMETERS']['ADVECTION SCHEME'] = 1
        self.cases.append(case_inputs)

        # HLL-SWE / SFZ's roots / Grass formula
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "bedload_soni.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'HLL-SVE'
        case_inputs['NUMERICAL PARAMETERS']['ADVECTION SCHEME'] = 2
        self.cases.append(case_inputs)

    def post(self, show=True):
        """
        Post
        """
        set_rcparams()
        c0, c1, c2 = get_default_color_palette()
        markers = get_default_markers()
        colors = c0 + c1 + c2
        fig, ax = plt.subplots(3, 2, figsize=(9., 8.))

        # Loop on cases:
        for i, case in enumerate(self.cases):
            RES = case['OUTPUT PARAMETERS']['RESULT FILE']
            file_ini = os.path.join(self.local_dir, 'RESU/'+RES+'ini.dat')
            file_fin = os.path.join(self.local_dir, 'RESU/'+RES+'fin.dat')
            ini = np.loadtxt(file_ini, delimiter=',', skiprows=2)
            res = np.loadtxt(file_fin, delimiter=',', skiprows=2)

            # extract values
            zb0= ini[:, 1]
            h0 = ini[:, 2]
            u0 = ini[:, 3]
            x = res[:, 0]
            zb= res[:, 1]
            h = res[:, 2]
            u = res[:, 3]
            z0 = res[:, 4]
            qs = res[:, 5]
            fr= abs(u)/np.sqrt(9.81*h)

            # plot h
            ax[0,0].plot(x, zb, color=colors[i], lw=0.75)
            #ax[0,0].fill_between(x, zb, min(zb), color='grey', alpha=0.5)

            ax[0,0].plot(x, zb+h, label=RES, color=colors[i], marker=markers[i], markersize=2, lw=1.)
            # plot q
            ax[0,1].plot(x, h*u, label=RES, color=colors[i], marker=markers[i], markersize=2, lw=1.)

            # plot b
            ax[1,0].plot(x, zb, label=RES, color=colors[i], marker=markers[i], markersize=2, lw=1.)
            # plot qs
            ax[1,1].plot(x, qs, label=RES, color=colors[i], marker=markers[i], markersize=2, lw=1.)            

            # plot Froude
            ax[2,0].plot(x, fr, label=RES, color=colors[i], marker=markers[i], markersize=2, lw=1.)            
            # plot u
            ax[2,1].plot(x, u, label=RES, color=colors[i], marker=markers[i], markersize=2, lw=1.)
        
        ax[0,0].legend(loc=2)
        ax[0,0].set_ylabel("$z_b+h$ (m)")
        ax[0,0].set_xlabel("$x$ (m)")
        ax[0,0].grid()
        ax[0,1].legend(loc=2)
        ax[0,1].set_ylabel("$q$ (m$^2$/s)")
        ax[0,1].set_xlabel("$x$ (m)")
        ax[0,1].grid()
        
        ax[1,1].legend(loc=2)
        ax[1,1].set_ylabel("$q_s$ (m$^2$/s)")
        ax[1,1].set_xlabel("$x$ (m)")
        ax[1,1].grid()
        ax[1,0].set_ylabel("$z_b$ (-)")
        ax[1,0].set_xlabel("$x$ (m)")
        ax[1,0].legend(loc=2)
        ax[1,0].grid()
        
        ax[2,1].legend(loc=2)
        ax[2,1].set_ylabel("$u$ (m/s)")
        ax[2,1].set_xlabel("$x$ (m)")
        ax[2,1].grid()
        ax[2,0].set_ylabel("$Fr$ (-)")
        ax[2,0].set_xlabel("$x$ (m)")
        ax[2,0].legend(loc=2)
        ax[2,0].grid()
        
        plt.savefig("FIG/{}.png".format(self.name), dpi=300)
        if show:
            plt.show()

        # load result
        res0 = np.loadtxt("RESU/ACU-SVEini.dat", delimiter=',', skiprows=2)
        x = res0[:, 0]
        zb0 = res0[:, 1]
        res = np.loadtxt("RESU/ACU-SVE_SNAP0001.dat", delimiter=',', skiprows=2)
        zb30 = res[:, 1]
        res = np.loadtxt("RESU/ACU-SVE_SNAP0002.dat", delimiter=',', skiprows=2)
        zb60 = res[:, 1]
        res = np.loadtxt("RESU/ACU-SVE_SNAP0003.dat", delimiter=',', skiprows=2)
        zb90 = res[:, 1]

        # load reference
        ref = np.loadtxt("REF/soni.txt")
        xr = ref[:, 0]
        zbr30 = ref[:, 1]
        zbr60 = ref[:, 2]
        zbr90 = ref[:, 3]

        # plot
        set_rcparams()
        c0, c1, c2 = get_default_color_palette()
        markers = get_default_markers()
        colors = c0 + c1 + c2
        fig, ax = plt.subplots(1, 1, figsize=(6., 5.))
        ax.plot(x, zb0, color='k', ls='--', lw=0.75, label='$z_b(0)$')
        ax.plot(x, zb30, color=c0[0], ls='-', lw=1.5, label='$z_b(t=30$min.$)$')
        ax.plot(x, zb60, color=c0[1], ls='-', lw=1.5, label='$z_b(t=60$min.$)$')
        ax.plot(x, zb90, color=c0[2], ls='-', lw=1.5, label='$z_b(t=90$min.$)$')
        ax.plot(xr, zbr30, color=c0[0], ls='-', marker='o', markersize=4, lw=0., label='Obs. $t=30$min.')
        ax.plot(xr, zbr60, color=c0[1], ls='-', marker='s', markersize=4, lw=0., label='Obs. $t=60$min.')
        ax.plot(xr, zbr90, color=c0[2], ls='-', marker='D', markersize=4, lw=0., label='Obs. $t=90$min.')
        ax.legend(loc=1)
        ax.set_ylim([-0.15, 0.05])
        ax.set_ylabel("$z$ (m)")
        ax.set_xlabel("$x$ (m)")
        ax.grid()
        plt.savefig("FIG/soni_MPM.png", dpi=300)
        if show:
            plt.show()


if __name__ == '__main__':

    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--run', default=False, action="store_true", help='run arrakis')
    parser.add_argument('-p', '--post', default=False, action="store_true", help='post processing')
    parser.add_argument('-c', '--check', default=False, action="store_true", help='check results')
    parser.add_argument('--reset-ref', default=False, action="store_true", help='WARNING: only use in extreme necessity')
    parser.add_argument('-nc', '--ncsize', type=int, default=4, help='number of OpenMP threads (default: 4)')
    args = parser.parse_args()

    # execute study
    study = bedload_soni_1()
    study.execute(args)
