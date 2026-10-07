#!/usr/bin/env python
import os
import argparse
import numpy as np
import matplotlib.pylab as plt
from arrakis import *

class TracerMultiple_1(Study):

    def pre(self):
        self.set_dir(file=__file__)
        self.name = "tracer_multiple"
        self.cases = []

        # HLL
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "tracer_multiple.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'HLL'
        case_inputs['NUMERICAL PARAMETERS']['SPACE ORDER'] = 1
        case_inputs['NUMERICAL PARAMETERS']['TIME SCHEME'] = 1
        case_inputs['NUMERICAL PARAMETERS']['SPACE ORDER FOR TRACERS'] = 1
        self.cases.append(case_inputs)

        # Steady kernel
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "tracer_multiple.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'STD'
        case_inputs['TIME PARAMETERS']['VARIABLE TIME STEP'] = False
        case_inputs['TIME PARAMETERS']['TIME STEP'] = 0.014
        case_inputs['NUMERICAL PARAMETERS']['TIME SCHEME'] = 0
        self.cases.append(case_inputs)


    def post(self, show=True):
        """
        Post
        """
        set_rcparams()
        c0, c1, c2 = get_default_color_palette()
        markers = get_default_markers()
        colors = c0 + c1 + c2
        fig, ax = plt.subplots(2, 1, figsize=(8, 5))

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
            T20 = restr10[:,5]
            T1 = restr1[:,4]
            T2 = restr1[:,5]
            
            # T ini
            if i==0:
                ax[0].plot(x, T10, label='T(t=0)', color='k', ls='--')
                ax[1].plot(x, T20, label='T(t=0)', color='k', ls='--')

            # plot
            ax[0].plot(x, T1, label=RES, color=c0[i], marker='o', markersize=2)
            ax[1].plot(x, T2, label=RES, color=c0[i], marker='s', markersize=2)

        ax[0].set_ylabel("$T_1$ (-)")
        ax[0].set_xlabel("$x$ (m)")
        ax[0].set_ylim([0.49, 0.51])
        ax[0].legend(loc=1)
        ax[0].grid()

        ax[1].set_ylabel("$T_2$ (-)")
        ax[1].set_xlabel("$x$ (m)")
        ax[1].set_ylim([0.99, 1.01])
        ax[1].legend(loc=1)
        ax[1].grid()
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
    study = TracerMultiple_1()
    study.execute(args)
