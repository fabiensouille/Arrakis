#!/usr/bin/env python
import os
import argparse
import numpy as np
import matplotlib.pylab as plt
from arrakis import *

class TracerCanal(Study):

    def pre(self):
        self.set_dir(file=__file__)
        self.name = "tracer_canal"
        self.cases = []

        # HLL
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "tracer_canal.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'HLL'
        case_inputs['TIME PARAMETERS']['VARIABLE TIME STEP'] = True
        case_inputs['NUMERICAL PARAMETERS']['TIME SCHEME'] = 1
        case_inputs['NUMERICAL PARAMETERS']['SPACE ORDER'] = 1
        case_inputs['NUMERICAL PARAMETERS']['SPACE ORDER FOR TRACERS'] = 1
        case_inputs['NUMERICAL PARAMETERS']['ADVECTION SCHEME'] = 2
        self.cases.append(case_inputs)

        # Steady kernel
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "tracer_canal.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'STD-CFL=1'
        case_inputs['TIME PARAMETERS']['VARIABLE TIME STEP'] = False
        case_inputs['TIME PARAMETERS']['TIME STEP'] = 0.014
        case_inputs['NUMERICAL PARAMETERS']['TIME SCHEME'] = 0
        self.cases.append(case_inputs)

        # Steady kernel
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "tracer_canal.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'STD-CFL=5'
        case_inputs['TIME PARAMETERS']['VARIABLE TIME STEP'] = False
        case_inputs['TIME PARAMETERS']['TIME STEP'] = 5*0.014
        case_inputs['NUMERICAL PARAMETERS']['TIME SCHEME'] = 0
        self.cases.append(case_inputs)

        # Steady kernel
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "tracer_canal.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'STD-CFL=10'
        case_inputs['TIME PARAMETERS']['VARIABLE TIME STEP'] = False
        case_inputs['TIME PARAMETERS']['TIME STEP'] = 10*0.014
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
        ax.legend()
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
    study = TracerCanal()
    study.execute(args)
