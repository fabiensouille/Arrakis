#!/usr/bin/env python
import os
import argparse
import numpy as np
import matplotlib.pylab as plt
from arrakis import *

class deposition(Study):

    def pre(self):
        self.set_dir(file=__file__)
        self.name = "deposition"
        self.cases = []

        # Ordre 1
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "deposition.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'DEP'
        case_inputs['NUMERICAL PARAMETERS']['TIME SCHEME'] = 1
        case_inputs['NUMERICAL PARAMETERS']['SPACE ORDER'] = 1
        case_inputs['NUMERICAL PARAMETERS']['SPACE ORDER FOR TRACERS'] = 1
        self.cases.append(case_inputs)
#        
#        # Ordre 2
#        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "deposition.yml"))
#        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'Ord.2'
#        case_inputs['NUMERICAL PARAMETERS']['TIME SCHEME'] = 2
#        case_inputs['NUMERICAL PARAMETERS']['SPACE ORDER'] = 2
#        case_inputs['NUMERICAL PARAMETERS']['SPACE ORDER FOR TRACERS'] = 2
#        self.cases.append(case_inputs)

#        # Ordre 5
#        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "deposition.yml"))
#        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'Ord.5'
#        case_inputs['NUMERICAL PARAMETERS']['TIME SCHEME'] = 2
#        case_inputs['NUMERICAL PARAMETERS']['SPACE ORDER'] = 2
#        case_inputs['NUMERICAL PARAMETERS']['SPACE ORDER FOR TRACERS'] = 5
#        self.cases.append(case_inputs)

    def post(self, show=True):
        """
        Post
        """
        set_rcparams()
        c0, c1, c2 = get_default_color_palette()
        markers = get_default_markers()
        colors = c0 + c1 + c2
        fig, ax = plt.subplots(2, 2, figsize=(9., 8.))

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
            zf = res[:, 4]
            qs = res[:, 5]
            C = res[:, 6]
            Ceq = res[:, 7]
            fr= abs(u)/np.sqrt(9.81*h)

            # plot h
            ax[0,0].plot(x, zb, color=colors[i], lw=0.75)
            ax[0,0].fill_between(x, zf, -0.6, color='saddlebrown', alpha=0.65)
            ax[0,0].fill_between(x, zb, zf, color='peru', alpha=0.45)
            ax[0,0].fill_between(x, zb, zb+h, color='steelblue', alpha=0.25)
            ax[0,0].plot(x, zb+h, label=RES, color=colors[i], marker=markers[i], markersize=2, lw=1.)
            
            # plot q
            ax[0,1].plot(x, h*u, label=RES, color=colors[i], marker=markers[i], markersize=2, lw=1.)

            # plot b
            ax[1,0].plot(x, zb, label=RES, color=colors[i], marker=markers[i], markersize=2, lw=1.)
            
            # plot C
            ax[1,1].plot(x, C, label=RES, color=colors[i], marker=markers[i], markersize=2, lw=1.)
            ax[1,1].plot(x, Ceq, label=RES+" (Ceq)", color='k', marker='o', markersize=2, lw=1.)

        ax[0,0].legend(loc=2)
        ax[0,0].set_ylabel("$z_b+h$ (m)")
        ax[0,0].set_xlabel("$x$ (m)")
        ax[0,0].grid()
        ax[0,1].legend(loc=2)
        ax[0,1].set_ylabel("$q$ (m$^2$/s)")
        ax[0,1].set_xlabel("$x$ (m)")
        ax[0,1].grid()
        
        ax[1,1].legend(loc=2)
        ax[1,1].set_ylabel("$C$ (kg/m$^3$)")
        ax[1,1].set_xlabel("$x$ (m)")
        ax[1,1].grid()
        ax[1,0].set_ylabel("$z_b$ (-)")
        ax[1,0].set_xlabel("$x$ (m)")
        ax[1,0].legend(loc=2)
        ax[1,0].grid()

        plt.savefig("FIG/{}.png".format(self.name), dpi=300)
        if show:
            plt.show()


    def post_trac(self, show=True):
        """
        Post
        """
        set_rcparams()
        c0, c1, c2 = get_default_color_palette()
        markers = get_default_markers()
        colors = c0 + c1 + c2
        fig, ax = plt.subplots(1, 1, figsize=(8, 3))

        # Loop on cases:
        for i, case in enumerate(self.cases):
            RES = case['OUTPUT PARAMETERS']['RESULT FILE']
            
            file_10 = os.path.join(self.local_dir, 'RESU/'+RES+'ini.dat')
            file_1 = os.path.join(self.local_dir, 'RESU/'+RES+'fin.dat')
            restr10 = np.loadtxt(file_10, delimiter=',', skiprows=2)
            restr1 = np.loadtxt(file_1, delimiter=',', skiprows=2)

            # arrays
            x = restr1[:,0]
            T10 = restr10[:,4]
            T1 = restr1[:,4]
            
            # T ini
            if i==0:
                ax.plot(x, T10, label='T(t=0)', color='k', ls='--')

            # plot
            ax.plot(x, T1, label=RES, color=c0[i], marker='o', markersize=2)

        ax.set_ylabel("$T_1$ (-)")
        ax.set_xlabel("$x$ (m)")
        ax.legend(loc=1)
        ax.grid()
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
    study = deposition()
    study.execute(args)
