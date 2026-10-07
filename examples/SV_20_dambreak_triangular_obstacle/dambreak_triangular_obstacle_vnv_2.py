#!/usr/bin/env python
import os
import argparse
import numpy as np
import matplotlib.pylab as plt
from arrakis import *

class dambreak_triangular_obstacle_2(Study):

    def pre(self):
        self.set_dir(file=__file__)
        self.name = "dambreak_triangular_obstacle_2"
        self.cases = []

        # Roe
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "dambreak_triangular_obstacle.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'ROE'
        case_inputs['OUTPUT PARAMETERS']['PROBES FILE NAME'] = 'PROBE-ROE'
        case_inputs['NUMERICAL PARAMETERS']['ADVECTION SCHEME'] = 1
        case_inputs['NUMERICAL PARAMETERS']['SPACE ORDER'] = 2
        case_inputs['NUMERICAL PARAMETERS']['TIME SCHEME'] = 2
        self.cases.append(case_inputs)

        # HLL
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "dambreak_triangular_obstacle.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'HLL'
        case_inputs['OUTPUT PARAMETERS']['PROBES FILE NAME'] = 'PROBE-HLL'
        case_inputs['NUMERICAL PARAMETERS']['ADVECTION SCHEME'] = 2
        case_inputs['NUMERICAL PARAMETERS']['SPACE ORDER'] = 2
        case_inputs['NUMERICAL PARAMETERS']['TIME SCHEME'] = 2
        self.cases.append(case_inputs)
        
        # Kinetic
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "dambreak_triangular_obstacle.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'KIN'
        case_inputs['OUTPUT PARAMETERS']['PROBES FILE NAME'] = 'PROBE-KIN'
        case_inputs['NUMERICAL PARAMETERS']['ADVECTION SCHEME'] = 4
        case_inputs['NUMERICAL PARAMETERS']['SPACE ORDER'] = 2
        case_inputs['NUMERICAL PARAMETERS']['TIME SCHEME'] = 2
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
            fr= res[:, 5]

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


        # plot probe results
        # ~~~~~~~~~~~~~~~~~~

        # load probe reference results
        r0 = np.loadtxt("REF/REF_G2.txt" , delimiter=',', skiprows=1)
        r1 = np.loadtxt("REF/REF_G4.txt" , delimiter=',', skiprows=1)
        r2 = np.loadtxt("REF/REF_G10.txt", delimiter=',', skiprows=1)
        r3 = np.loadtxt("REF/REF_G13.txt", delimiter=',', skiprows=1)
        r4 = np.loadtxt("REF/REF_G20.txt", delimiter=',', skiprows=1)

        # load probe results
        p0 = np.loadtxt("RESU/PROBE-HLL_0001.dat", delimiter=',', skiprows=2)
        p1 = np.loadtxt("RESU/PROBE-HLL_0002.dat", delimiter=',', skiprows=2)
        p2 = np.loadtxt("RESU/PROBE-HLL_0003.dat", delimiter=',', skiprows=2)
        p3 = np.loadtxt("RESU/PROBE-HLL_0004.dat", delimiter=',', skiprows=2)
        p4 = np.loadtxt("RESU/PROBE-HLL_0005.dat", delimiter=',', skiprows=2)

        # plot p1
        set_rcparams()
        c0, c1, c2 = get_default_color_palette()
        markers = get_default_markers()
        colors = c0 + c1 + c2
        fig, ax = plt.subplots(1, 1, figsize=(5., 4.))
        ax.plot(r1[:,0], r1[:,1], label='Obs.', color=c0[1], lw=0., marker='o', markersize=4)
        ax.plot(p1[:,0], p1[:,1], label='Simulation', color=c0[0], lw=2.)
        ax.set_ylabel("$h$ (m)")
        ax.set_xlabel("$t$ (s)")
        plt.text(2.5, 0.52, '$G4$', fontdict=None, )
        plt.xlim([0., 40.])
        plt.ylim([0., 0.6])
        plt.legend()
        plt.grid()
        plt.savefig("FIG/vnv2_probe_1.png", dpi=300)
        if show:
            plt.show()
        
        # plot p2
        set_rcparams()
        c0, c1, c2 = get_default_color_palette()
        markers = get_default_markers()
        colors = c0 + c1 + c2
        fig, ax = plt.subplots(1, 1, figsize=(5., 4.))
        ax.plot(r2[:,0], r2[:,1], label='Obs.', color=c0[1], lw=0., marker='o', markersize=4)
        ax.plot(p2[:,0], p2[:,1], label='Simulation', color=c0[0], lw=2.)
        ax.set_ylabel("$h$ (m)")
        ax.set_xlabel("$t$ (s)")
        plt.text(2.5, 0.72, '$G10$', fontdict=None, )
        plt.xlim([0., 40.])
        plt.ylim([0., 0.8])
        plt.legend()
        plt.grid()
        plt.savefig("FIG/vnv2_probe_2.png", dpi=300)
        if show:
            plt.show()
        
        # plot p3
        set_rcparams()
        c0, c1, c2 = get_default_color_palette()
        markers = get_default_markers()
        colors = c0 + c1 + c2
        fig, ax = plt.subplots(1, 1, figsize=(5., 4.))
        ax.plot(r3[:,0], r3[:,1], label='Obs.', color=c0[1], lw=0., marker='o', markersize=4)
        ax.plot(p3[:,0], p3[:,1], label='Simulation', color=c0[0], lw=2.)
        ax.set_ylabel("$h$ (m)")
        ax.set_xlabel("$t$ (s)")
        plt.text(2.5, 0.22, '$G13$', fontdict=None, )
        plt.xlim([0., 40.])
        plt.ylim([0., 0.3])
        plt.legend()
        plt.grid()
        plt.savefig("FIG/vnv2_probe_3.png", dpi=300)
        if show:
            plt.show()
        
        # plot p4
        set_rcparams()
        c0, c1, c2 = get_default_color_palette()
        markers = get_default_markers()
        colors = c0 + c1 + c2
        fig, ax = plt.subplots(1, 1, figsize=(5., 4.))
        ax.plot(r4[:,0], r4[:,1], label='Obs.', color=c0[1], lw=0., marker='o', markersize=4)
        ax.plot(p4[:,0], p4[:,1], label='Simulation', color=c0[0], lw=2.)
        ax.set_ylabel("$h$ (m)")
        ax.set_xlabel("$t$ (s)")
        plt.text(2.5, 0.52, '$G20$', fontdict=None, )
        plt.xlim([0., 40.])
        plt.ylim([0., 0.6])
        plt.legend()
        plt.grid()
        plt.savefig("FIG/vnv2_probe_4.png", dpi=300)
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
    study = dambreak_triangular_obstacle_2()
    study.execute(args)
