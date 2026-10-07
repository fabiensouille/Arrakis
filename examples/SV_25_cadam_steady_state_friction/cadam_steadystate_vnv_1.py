#!/usr/bin/env python
import os
import argparse
import numpy as np
import matplotlib.pylab as plt
from arrakis import *

class steadystate_1(Study):

    def pre(self):
        self.set_dir(file=__file__)
        self.name = "steadystate_1"
        self.cases = []

        # HLL + HR
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "cadam_steadystate_1.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'HLL-HR[zb]-CEN[B]'
        case_inputs['NUMERICAL PARAMETERS']['ADVECTION SCHEME'] = 2
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['UNLOCK ADVANCED PARAMETERS'] = True
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['HYDROSTATIC RECONSTRUCTION'] = 1
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['VARIABLE WIDTH SCHEME'] = 2
        self.cases.append(case_inputs)

        # Roe
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "cadam_steadystate_1.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'ROE-HR[zb]-CEN[B]'
        case_inputs['NUMERICAL PARAMETERS']['ADVECTION SCHEME'] = 1
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['UNLOCK ADVANCED PARAMETERS'] = True
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['HYDROSTATIC RECONSTRUCTION'] = 1
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['VARIABLE WIDTH SCHEME'] = 2
        self.cases.append(case_inputs)
        
        # Roe
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "cadam_steadystate_1.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'ROE-VC[zb]-CEN[B]'
        case_inputs['NUMERICAL PARAMETERS']['ADVECTION SCHEME'] = 1
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['UNLOCK ADVANCED PARAMETERS'] = True
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['HYDROSTATIC RECONSTRUCTION'] = 5
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['VARIABLE WIDTH SCHEME'] = 2
        self.cases.append(case_inputs)

        # Roe
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "cadam_steadystate_1.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'ROE-VC[zb]-VC[B]'
        case_inputs['NUMERICAL PARAMETERS']['ADVECTION SCHEME'] = 1
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['UNLOCK ADVANCED PARAMETERS'] = True
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['HYDROSTATIC RECONSTRUCTION'] = 5
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['VARIABLE WIDTH SCHEME'] = 5
        self.cases.append(case_inputs)
        
        # Roe
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "cadam_steadystate_1.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'ROE-HR[zb]-VC[B]'
        case_inputs['NUMERICAL PARAMETERS']['ADVECTION SCHEME'] = 1
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['UNLOCK ADVANCED PARAMETERS'] = True
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['HYDROSTATIC RECONSTRUCTION'] = 1
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['VARIABLE WIDTH SCHEME'] = 5
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
            q = res[:, 4]
            fr = np.zeros(len(u), dtype='d')
            for j in range(len(u)):
                fr[j] = abs(u[j])/np.sqrt(9.81*max(1.e-16, h[j]))

            # plot h
            if i==0:
                ax[0,0].plot(x, zb, label='z', color='k', lw=0.5)
                ax[0,0].fill_between(x, zb, min(zb), color='grey', alpha=0.5)
            ax[0,0].plot(x, zb+h, label=RES, color=colors[i], marker=markers[i], markersize=2, lw=1.)
            # plot q
            ax[0,1].plot(x, q, label=RES, color=colors[i], marker=markers[i], markersize=2, lw=1.)
            # plot u
            ax[1,1].plot(x, u, label=RES, color=colors[i], marker=markers[i], markersize=2, lw=1.)
            # plot Froude
            ax[1,0].plot(x, fr, label=RES, color=colors[i], marker=markers[i], markersize=2, lw=1.)

        # analytical solution
        ref_data = np.loadtxt('INI_149.ini', skiprows=2)
        ref_x = ref_data[:, 0]
        ref_h = ref_data[:, 1]
        ref_u = ref_data[:, 2]
        ax[0,0].plot(ref_x, zb+ref_h, label="Exact", color="k", ls='--', lw=1.)
        ax[0,1].plot(ref_x, ref_u*ref_h, label="Exact", color="k", ls='--', lw=1.)
        ax[1,1].plot(ref_x, ref_u, label="Exact", color="k", ls='--', lw=1.)

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

    # execute study
    study = steadystate_1()
    study.execute(args)
