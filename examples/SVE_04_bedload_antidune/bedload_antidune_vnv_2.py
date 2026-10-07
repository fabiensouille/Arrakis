#!/usr/bin/env python
import os
import argparse
import numpy as np
import matplotlib.pylab as plt
from arrakis import *

class bedload_antidune_2(Study):

    def pre(self):
        self.set_dir(file=__file__)
        self.name = "bedload_antidune_2"
        self.cases = []

        # ROE-SWE / Cardano's roots / Grass formula
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "bedload_antidune.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'ROE-SVE'
        case_inputs['MESH PARAMETERS']['NUMBER OF NODES'] = 201
        case_inputs['NUMERICAL PARAMETERS']['ADVECTION SCHEME'] = 1
        case_inputs['NUMERICAL PARAMETERS']['SPACE ORDER'] = 2
        case_inputs['NUMERICAL PARAMETERS']['TIME SCHEME'] = 2
        case_inputs['SEDIMENT PARAMETERS']['BEDLOAD FORMULA'] = 1
        self.cases.append(case_inputs)

        # HLL-SWE / SFZ's roots / Grass formula
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "bedload_antidune.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'HLL-SVE'
        case_inputs['MESH PARAMETERS']['NUMBER OF NODES'] = 201
        case_inputs['NUMERICAL PARAMETERS']['ADVECTION SCHEME'] = 2
        case_inputs['NUMERICAL PARAMETERS']['SPACE ORDER'] = 2
        case_inputs['NUMERICAL PARAMETERS']['TIME SCHEME'] = 2
        case_inputs['SEDIMENT PARAMETERS']['BEDLOAD FORMULA'] = 1
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
    study = bedload_antidune_2()
    study.execute(args)
