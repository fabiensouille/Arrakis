#!/usr/bin/env python
import os
import argparse
import numpy as np
import matplotlib.pylab as plt
from arrakis import *

class vanrijn(Study):

    def pre(self):
        self.set_dir(file=__file__)
        self.name = "vanrijn_calibration_VR"
        self.cases = []

        # Final Time
        self.tf = 18000

        # Friction
        #  3: Strickler (1981): ks/h in [0.007;0.1]
        #  4: Nikuradse (1950): ks/h in [0;0.1]
        #  5: Colebrook-White (1939): ks/h in [0;0.1]
        #  6: Haaland (1983) : explicit version of CW (C=12)
        #  7: Barr (1977)    : explicit version of CW (C=14.8)
        #  8: Bathurst (1985): ks/h in [0.1;5.15]
        #  9: Machiels (2011): Barr+Bathurst, ks/h in [0,5.15]
        friction = 4 
        friction_coef = 3.2e-2 # eq. Strickler = 46 (Ramette)
        
        # toggle rouse profile
        use_rouse_profile = True

        # Van Rijn (1984)
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "vanrijn.yml"))
        case_inputs['TIME PARAMETERS']['DURATION'] = self.tf
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'VR-DEF'
        case_inputs['PHYSICAL PARAMETERS']['BOTTOM FRICTION LAW'] = friction
        case_inputs['PHYSICAL PARAMETERS']['BOTTOM FRICTION COEFFICIENT'] = friction_coef
        case_inputs['SEDIMENT PARAMETERS']['SUSPENSION FORMULA'] = 3
        case_inputs['SEDIMENT PARAMETERS']['SUSPENSION CALIBRATION COEFFICIENT'] = 1.
        case_inputs['SEDIMENT PARAMETERS']['SUSPENSION SETTLING VELOCITY'] = 0.065
        case_inputs['SEDIMENT PARAMETERS']['SUSPENSION USE ROUSE PROFILE'] = use_rouse_profile
        case_inputs['SEDIMENT PARAMETERS']['BED FRICTION SCALING COEFFICIENT'] = 1.
        self.cases.append(case_inputs)

        # Van Rijn (1984)
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "vanrijn.yml"))
        case_inputs['TIME PARAMETERS']['DURATION'] = self.tf
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'VR-ws*1p5-ks*0p8'
        case_inputs['PHYSICAL PARAMETERS']['BOTTOM FRICTION LAW'] = friction
        case_inputs['PHYSICAL PARAMETERS']['BOTTOM FRICTION COEFFICIENT'] = friction_coef
        case_inputs['SEDIMENT PARAMETERS']['SUSPENSION FORMULA'] = 3
        case_inputs['SEDIMENT PARAMETERS']['SUSPENSION CALIBRATION COEFFICIENT'] = 1.5
        case_inputs['SEDIMENT PARAMETERS']['SUSPENSION SETTLING VELOCITY'] = 0.065
        case_inputs['SEDIMENT PARAMETERS']['SUSPENSION USE ROUSE PROFILE'] = use_rouse_profile
        case_inputs['SEDIMENT PARAMETERS']['BED FRICTION SCALING COEFFICIENT'] = 0.8
        self.cases.append(case_inputs)

        # Van Rijn (1984)
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "vanrijn.yml"))
        case_inputs['TIME PARAMETERS']['DURATION'] = self.tf
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'VR-ws*2p0-ks*0p6'
        case_inputs['PHYSICAL PARAMETERS']['BOTTOM FRICTION LAW'] = friction
        case_inputs['PHYSICAL PARAMETERS']['BOTTOM FRICTION COEFFICIENT'] = friction_coef
        case_inputs['SEDIMENT PARAMETERS']['SUSPENSION FORMULA'] = 3
        case_inputs['SEDIMENT PARAMETERS']['SUSPENSION CALIBRATION COEFFICIENT'] = 2.
        case_inputs['SEDIMENT PARAMETERS']['SUSPENSION SETTLING VELOCITY'] = 0.065
        case_inputs['SEDIMENT PARAMETERS']['SUSPENSION USE ROUSE PROFILE'] = use_rouse_profile
        case_inputs['SEDIMENT PARAMETERS']['BED FRICTION SCALING COEFFICIENT'] = 0.6
        self.cases.append(case_inputs)

        # Van Rijn (1984)
        case_inputs = read_yml_input_file(os.path.join(self.local_dir, "vanrijn.yml"))
        case_inputs['TIME PARAMETERS']['DURATION'] = self.tf
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'VR-ws*2p5-ks*0p4'
        case_inputs['PHYSICAL PARAMETERS']['BOTTOM FRICTION LAW'] = friction
        case_inputs['PHYSICAL PARAMETERS']['BOTTOM FRICTION COEFFICIENT'] = friction_coef
        case_inputs['SEDIMENT PARAMETERS']['SUSPENSION FORMULA'] = 3
        case_inputs['SEDIMENT PARAMETERS']['SUSPENSION CALIBRATION COEFFICIENT'] = 2.5
        case_inputs['SEDIMENT PARAMETERS']['SUSPENSION SETTLING VELOCITY'] = 0.065
        case_inputs['SEDIMENT PARAMETERS']['SUSPENSION USE ROUSE PROFILE'] = use_rouse_profile
        case_inputs['SEDIMENT PARAMETERS']['BED FRICTION SCALING COEFFICIENT'] = 0.4
        self.cases.append(case_inputs)


    def post(self, show=True):
        """
        Post
        """
        set_rcparams()
        c0, c1, c2 = get_default_color_palette()
        markers = get_default_markers()
        colors = c0 + c1 + c2
        fig, ax = plt.subplots(2, 2, figsize=(16., 10.))

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
            ax[0,0].fill_between(x, zf, -0.4, color='saddlebrown', alpha=0.65)
            ax[0,0].fill_between(x, zb, zf, color='peru', alpha=0.45)
            ax[0,0].fill_between(x, zb, zb+h, color='steelblue', alpha=0.25)
            ax[0,0].plot(x, zb+h, label=RES, color=colors[i], marker=markers[i], markersize=2, lw=1.)
            
            # plot q
            ax[0,1].plot(x, h*u, label=RES, color=colors[i], marker=markers[i], markersize=2, lw=1.)

            # plot b
            ax[1,0].plot(x, zb, label=RES, color=colors[i], marker=markers[i], markersize=2, lw=1.)
            
            # plot C
            ax[1,1].plot(x, C, label=RES, color=colors[i], marker=markers[i], markersize=0, lw=1.)
            ax[1,1].plot(x, Ceq, label=RES+" (Ceq)", color=colors[i], ls='--', marker='o', markersize=2, lw=1., markevery=10)

        # reference
        ref = np.loadtxt("REF/OBS.txt", delimiter=',')
        ax[1,0].plot(x, zb0, c='k', ls='--', lw=1, label="$z_b(0)$ (m)") 
        ax[1,0].plot(ref[:,0], ref[:,1], c='k', marker='o', markersize=4, lw=0, label="$Obs.$") 

        ax[0,0].legend(loc=2)
        ax[0,0].set_ylabel("$z_b+h$ (m)")
        ax[0,0].set_xlabel("$x$ (m)")
        ax[0,0].grid()

        ax[0,1].legend(loc=2)
        ax[0,1].set_ylabel("$q$ (m$^2$/s)")
        ax[0,1].set_xlabel("$x$ (m)")
        ax[0,1].grid()
        
        ax[1,1].legend(loc=2)
        ax[1,1].set_ylabel("$C$ ($g/L$)")
        ax[1,1].set_xlabel("$x$ (m)")
        ax[1,1].grid()

        ax[1,0].set_ylabel("$z_b$ (-)")
        ax[1,0].set_xlabel("$x$ (m)")
        ax[1,0].legend(loc=2)
        ax[1,0].set_ylim([-0.2, 0.1])
        ax[1,0].grid()

        plt.suptitle("Result at time = ${} s$".format(self.tf))
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
    study = vanrijn()
    study.execute(args)
