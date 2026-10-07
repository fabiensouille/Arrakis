#!/usr/bin/env python
import os
import argparse
import numpy as np
import matplotlib.pylab as plt
from arrakis import *

class bumptrans_1(Study):

    def pre(self):
        self.set_dir(file=__file__)
        self.name = "bumptrans_1"
        self.cases = []

        # Roe
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "bumptrans.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'ROE'
        case_inputs['NUMERICAL PARAMETERS']['ADVECTION SCHEME'] = 1
        case_inputs['NUMERICAL PARAMETERS']['HYDROSTATIC RECONSTRUCTION'] = 1
        self.cases.append(case_inputs)

        # HLL
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "bumptrans.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'HLL'
        case_inputs['NUMERICAL PARAMETERS']['ADVECTION SCHEME'] = 2
        case_inputs['NUMERICAL PARAMETERS']['HYDROSTATIC RECONSTRUCTION'] = 1
        self.cases.append(case_inputs)

        # Kinetic
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "bumptrans.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'KIN'
        case_inputs['NUMERICAL PARAMETERS']['ADVECTION SCHEME'] = 4
        case_inputs['NUMERICAL PARAMETERS']['HYDROSTATIC RECONSTRUCTION'] = 1
        self.cases.append(case_inputs)


    def post(self, show=True):
        """
        Post
        """
        set_rcparams()
        c0, c1, c2 = get_default_color_palette()
        markers = get_default_markers()
        colors = c0 + c1 + c2
        fig, ax = plt.subplots(2, 2, figsize=(8.5, 8.))

        # analytical solution
        sol = BumpAnalyticSol(flow='trans', Q=0.4, hl=0.4, bottom_function='exponential')

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
            fr= abs(u)/np.sqrt(9.81*h)

            # plot h
            if i==0:
                ax[0,0].plot(x, zb, label='z', color='k', lw=0.5)
                ax[0,0].fill_between(x, zb, min(zb), color='grey', alpha=0.5)
            ax[0,0].plot(x, zb+h, label=RES, color=colors[i], marker=markers[i], markersize=2, lw=1.)
            # plot q
            ax[0,1].plot(x, h*u, label=RES, color=colors[i], marker=markers[i], markersize=2, lw=1.)
            # plot u
            ax[1,1].plot(x, u, label=RES, color=colors[i], marker=markers[i], markersize=2, lw=1.)
            # plot Froude
            ax[1,0].plot(x, fr, label=RES, color=colors[i], marker=markers[i], markersize=2, lw=1.)
        
        # analytical solution
        sol()
        ref = [sol.x, sol.zb, sol.H, sol.U]
        ax[0,0].plot(ref[0], ref[1]+ref[2], label="Exact", color="k", ls='--', lw=1.)
        ax[0,1].plot(ref[0], ref[2]*ref[3], label="Exact", color="k", ls='--', lw=1.)
        
        ax[0,0].legend(loc=2)
        ax[0,0].set_ylabel("$z$ (m)")
        ax[0,0].set_xlabel("$x$ (m)")
        ax[0,0].grid()
        ax[0,1].legend(loc=2)
        ax[0,1].set_ylabel("$q$ (m$^2$/s)")
        ax[0,1].set_xlabel("$x$ (m)")
        ax[0,1].grid()
        ax[1,1].legend(loc=2)
        ax[1,1].set_ylabel("$u$ (m/s)")
        ax[1,1].set_xlabel("$x$ (m)")
        ax[1,1].grid()
        ax[1,0].set_ylabel("$Fr$ (-)")
        ax[1,0].set_xlabel("$x$ (m)")
        ax[1,0].legend(loc=1)
        ax[1,0].grid()
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
    
    if args.run==False and args.post==False and args.check==False and args.reset_ref==False:
        no_args = True
    else:
        no_args = False

    # execute study
    study = bumptrans_1()
    study.execute(args)
