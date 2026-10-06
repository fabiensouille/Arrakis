#!/usr/bin/env python
import os
import time
import numpy as np
import matplotlib.pylab as plt
from arrakis import *

ROOT_DIR = os.path.dirname(__file__)

class VNV_bedload_dambreak_qs():

    def vnv_pre(self):
        self.name = "VNV_bedload_dambreak_qs"
        self.cases = []

        # GLOBAL PARAMETERS
        Ks = 35.
        rhos = 2600.
        d50 = 1.e-3
        mpm_tausc = 0.047
        
        # analogy with MPM to get Ag
        g = 9.8065
        Rh = 0.2
        R = (rhos/1000.)-1.
        grass_ag = 8.*np.sqrt(g)/(R*(Ks**3.)*np.sqrt(Rh))
        grass_mg = 3.
        print("Ag =", grass_ag)

        # GRASS
        # *****
        # Grass formula / Exact derivatives
        case_inputs = read_yml_input_file(os.path.join(ROOT_DIR, "bedload_dambreak_qs.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'Grass-exact'
        case_inputs['PHYSICAL PARAMETERS']['BOTTOM FRICTION COEFFICIENT'] = Ks
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['BEDLOAD FLUX DERIVATIVES'] = 0
        case_inputs['SEDIMENT PARAMETERS']['BEDLOAD FORMULA'] = 1
        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT DENSITY'] = rhos
        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT D50'] = d50
        case_inputs['SEDIMENT PARAMETERS']['GRASS AG'] = grass_ag
        case_inputs['SEDIMENT PARAMETERS']['GRASS MG'] = grass_mg
        self.cases.append(case_inputs)
        
        # Grass formula / Approx derivatives
        case_inputs = read_yml_input_file(os.path.join(ROOT_DIR, "bedload_dambreak_qs.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'Grass-approx'
        case_inputs['PHYSICAL PARAMETERS']['BOTTOM FRICTION COEFFICIENT'] = Ks
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['BEDLOAD FLUX DERIVATIVES'] = 4
        case_inputs['SEDIMENT PARAMETERS']['BEDLOAD FORMULA'] = 1
        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT DENSITY'] = rhos
        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT D50'] = d50
        case_inputs['SEDIMENT PARAMETERS']['GRASS AG'] = grass_ag
        case_inputs['SEDIMENT PARAMETERS']['GRASS MG'] = grass_mg
        self.cases.append(case_inputs)

        # MPM
        # ***
        # MPM formula / Exact derivatives
        case_inputs = read_yml_input_file(os.path.join(ROOT_DIR, "bedload_dambreak_qs.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'MPM-exact'
        case_inputs['PHYSICAL PARAMETERS']['BOTTOM FRICTION COEFFICIENT'] = Ks
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['BEDLOAD FLUX DERIVATIVES'] = 0
        case_inputs['SEDIMENT PARAMETERS']['BEDLOAD FORMULA'] = 3
        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT DENSITY'] = rhos
        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT D50'] = d50
        case_inputs['SEDIMENT PARAMETERS']['CRITICAL SHIELDS'] = mpm_tausc
        self.cases.append(case_inputs)
        
        # MPM formula / Approx derivatives
        case_inputs = read_yml_input_file(os.path.join(ROOT_DIR, "bedload_dambreak_qs.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'MPM-approx'
        case_inputs['PHYSICAL PARAMETERS']['BOTTOM FRICTION COEFFICIENT'] = Ks
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['BEDLOAD FLUX DERIVATIVES'] = 4
        case_inputs['SEDIMENT PARAMETERS']['BEDLOAD FORMULA'] = 3
        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT DENSITY'] = rhos
        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT D50'] = d50
        case_inputs['SEDIMENT PARAMETERS']['CRITICAL SHIELDS'] = mpm_tausc
        self.cases.append(case_inputs)
        
        # STANDARD MPM
        # ************
        # standardized MPM formula / Exact derivatives
        case_inputs = read_yml_input_file(os.path.join(ROOT_DIR, "bedload_dambreak_qs.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'sMPM-exact'
        case_inputs['PHYSICAL PARAMETERS']['BOTTOM FRICTION COEFFICIENT'] = Ks
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['BEDLOAD FLUX DERIVATIVES'] = 0
        case_inputs['SEDIMENT PARAMETERS']['BEDLOAD FORMULA'] = 4
        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT DENSITY'] = rhos
        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT D50'] = d50
        case_inputs['SEDIMENT PARAMETERS']['CRITICAL SHIELDS'] = mpm_tausc
        case_inputs['SEDIMENT PARAMETERS']['STD MPM COEF A'] = 10.24
        case_inputs['SEDIMENT PARAMETERS']['STD MPM COEF B'] = 0.7
        self.cases.append(case_inputs)
        
        # standardized MPM formula / Approx derivatives
        case_inputs = read_yml_input_file(os.path.join(ROOT_DIR, "bedload_dambreak_qs.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'sMPM-approx'
        case_inputs['PHYSICAL PARAMETERS']['BOTTOM FRICTION COEFFICIENT'] = Ks
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['BEDLOAD FLUX DERIVATIVES'] = 4
        case_inputs['SEDIMENT PARAMETERS']['BEDLOAD FORMULA'] = 4
        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT DENSITY'] = rhos
        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT D50'] = d50
        case_inputs['SEDIMENT PARAMETERS']['CRITICAL SHIELDS'] = mpm_tausc
        case_inputs['SEDIMENT PARAMETERS']['STD MPM COEF A'] = 10.24
        case_inputs['SEDIMENT PARAMETERS']['STD MPM COEF B'] = 0.7
        self.cases.append(case_inputs)

        # Van Rijn
        # ********
        # Engelund formula / Exact derivatives
        case_inputs = read_yml_input_file(os.path.join(ROOT_DIR, "bedload_dambreak_qs.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'VR-exact'
        case_inputs['PHYSICAL PARAMETERS']['BOTTOM FRICTION COEFFICIENT'] = Ks
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['BEDLOAD FLUX DERIVATIVES'] = 0
        case_inputs['SEDIMENT PARAMETERS']['BEDLOAD FORMULA'] = 9
        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT DENSITY'] = rhos
        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT D50'] = d50
        self.cases.append(case_inputs)

        # Van Rijn
        # ********
        # Engelund formula / Approx derivatives
        case_inputs = read_yml_input_file(os.path.join(ROOT_DIR, "bedload_dambreak_qs.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'VR-approx'
        case_inputs['PHYSICAL PARAMETERS']['BOTTOM FRICTION COEFFICIENT'] = Ks
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['BEDLOAD FLUX DERIVATIVES'] = 1
        case_inputs['SEDIMENT PARAMETERS']['BEDLOAD FORMULA'] = 9
        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT DENSITY'] = rhos
        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT D50'] = d50
        self.cases.append(case_inputs)

#        # Engelund
#        # ********
#        # Engelund formula / Exact derivatives
#        case_inputs = read_yml_input_file(os.path.join(ROOT_DIR, "bedload_dambreak_qs.yml"))
#        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'EH-exact'
#        case_inputs['PHYSICAL PARAMETERS']['BOTTOM FRICTION COEFFICIENT'] = Ks
#        case_inputs['ADVANCED NUMERICAL PARAMETERS']['BEDLOAD FLUX DERIVATIVES'] = 0
#        case_inputs['SEDIMENT PARAMETERS']['BEDLOAD FORMULA'] = 8
#        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT DENSITY'] = rhos
#        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT D50'] = d50
#        self.cases.append(case_inputs)

#        # Engelund
#        # ********
#        # Engelund formula / Approx derivatives
#        case_inputs = read_yml_input_file(os.path.join(ROOT_DIR, "bedload_dambreak_qs.yml"))
#        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'EH-approx'
#        case_inputs['PHYSICAL PARAMETERS']['BOTTOM FRICTION COEFFICIENT'] = Ks
#        case_inputs['ADVANCED NUMERICAL PARAMETERS']['BEDLOAD FLUX DERIVATIVES'] = 1
#        case_inputs['SEDIMENT PARAMETERS']['BEDLOAD FORMULA'] = 8
#        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT DENSITY'] = rhos
#        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT D50'] = d50
#        self.cases.append(case_inputs)

#        # Lefort
#        # ******
#        # Lefort formula / Approx derivatives
#        case_inputs = read_yml_input_file(os.path.join(ROOT_DIR, "bedload_dambreak_qs.yml"))
#        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'Lefort-approx'
#        case_inputs['PHYSICAL PARAMETERS']['BOTTOM FRICTION COEFFICIENT'] = Ks
#        case_inputs['ADVANCED NUMERICAL PARAMETERS']['BEDLOAD FLUX DERIVATIVES'] = 1
#        case_inputs['SEDIMENT PARAMETERS']['BEDLOAD FORMULA'] = 5
#        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT DENSITY'] = rhos
#        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT D50'] = d50
#        self.cases.append(case_inputs)

#        # Recking 2013
#        # ************
#        # Lefort formula / Approx derivatives
#        case_inputs = read_yml_input_file(os.path.join(ROOT_DIR, "bedload_dambreak_qs.yml"))
#        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'Recking2013-approx'
#        case_inputs['PHYSICAL PARAMETERS']['BOTTOM FRICTION COEFFICIENT'] = Ks
#        case_inputs['ADVANCED NUMERICAL PARAMETERS']['BEDLOAD FLUX DERIVATIVES'] = 1
#        case_inputs['SEDIMENT PARAMETERS']['BEDLOAD FORMULA'] = 6
#        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT DENSITY'] = rhos
#        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT D50'] = d50
#        self.cases.append(case_inputs)

#        # Recking 2015
#        # ************
#        # Lefort formula / Approx derivatives
#        case_inputs = read_yml_input_file(os.path.join(ROOT_DIR, "bedload_dambreak_qs.yml"))
#        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'Recking2015-approx'
#        case_inputs['PHYSICAL PARAMETERS']['BOTTOM FRICTION COEFFICIENT'] = Ks
#        case_inputs['ADVANCED NUMERICAL PARAMETERS']['BEDLOAD FLUX DERIVATIVES'] = 1
#        case_inputs['SEDIMENT PARAMETERS']['BEDLOAD FORMULA'] = 7
#        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT DENSITY'] = rhos
#        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT D50'] = d50
#        self.cases.append(case_inputs)

    def vnv_run(self):
        """
        Run 
        """
        rep = os.path.join(ROOT_DIR, "RESU")
        os.system("rm {}/*.dat".format(rep))
        for case in self.cases:
            t0 = time.time()
            run_arrakis(case)
            t1 = time.time()
            print("CPU Time:", t1-t0 )

    def vnv_post_debug(self, show=True):
        """
        Post
        """
        set_rcparams()
        c0, c1, c2 = get_default_color_palette()
        markers = get_default_markers()
        markers += markers
        colors = c0 + c1 + c2 + ['m', 'c']
        fig, ax = plt.subplots(3, 2, figsize=(10., 9.))
        lsl = ['-', '--','-', '--','-', '--','--', '--','--', '--','--', '--']

        # Loop on cases:
        for i, case in enumerate(self.cases):
            RES = case['OUTPUT PARAMETERS']['RESULT FILE']
            file_ini = os.path.join(ROOT_DIR, 'RESU/'+RES+'ini.dat')
            file_fin = os.path.join(ROOT_DIR, 'RESU/'+RES+'fin.dat')
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
            taus = res[:, 6]
            fr= abs(u)/np.sqrt(9.81*h)

            # plot h
            ax[0,0].plot(x, zb, color=colors[i], lw=0.75)
            ax[0,0].plot(x, zb+h, label=RES, color=colors[i], marker=markers[i], 
                         markersize=3, lw=1., ls=lsl[i], markevery=8)
            # plot q
            ax[1,0].plot(x, h*u, label=RES, color=colors[i], marker=markers[i],
                         markersize=3, lw=1., ls=lsl[i], markevery=8)
            # plot b
            ax[0,1].plot(x, zb, label=RES, color=colors[i], marker=markers[i], 
                         markersize=3, lw=1., ls=lsl[i], markevery=8)
            # plot qs
            ax[1,1].plot(x, qs, label=RES, color=colors[i], marker=markers[i], 
                         markersize=3, lw=1., ls=lsl[i], markevery=8)
            # plot Froude
            ax[2,0].plot(x, fr, label=RES, color=colors[i], marker=markers[i], 
                         markersize=3, lw=1., ls=lsl[i], markevery=8)
            # plot taus
            ax[2,1].plot(x, taus, label=RES, color=colors[i], marker=markers[i], 
                         markersize=3, lw=1., ls=lsl[i], markevery=8)
        
        ax[0,0].legend(loc=2)
        ax[0,0].set_ylabel("$z_b+h$ (m)")
        ax[0,0].set_xlabel("$x$ (m)")
        ax[0,0].grid()
        
        ax[1,0].legend(loc=2)
        ax[1,0].set_ylabel("$q$ (m$^2$/s)")
        ax[1,0].set_xlabel("$x$ (m)")
        ax[1,0].grid()
        
        ax[1,1].legend(loc=2)
        ax[1,1].set_ylabel("$q_s$ (m$^2$/s)")
        ax[1,1].set_xlabel("$x$ (m)")
        ax[1,1].grid()
        
        ax[0,1].set_ylabel("$z_b$ (-)")
        ax[0,1].set_xlabel("$x$ (m)")
        ax[0,1].legend(loc=2)
        ax[0,1].grid()
        
        ax[2,1].legend(loc=2)
        ax[2,1].set_ylabel("$\\tau^*$")
        ax[2,1].set_xlabel("$x$ (m)")
        ax[2,1].grid()
        
        ax[2,0].set_ylabel("$Fr$ (-)")
        ax[2,0].set_xlabel("$x$ (m)")
        ax[2,0].legend(loc=2)
        ax[2,0].grid()
        
        plt.savefig("FIG/bedload_dambreak_qs_debug.png", dpi=300)
        if show:
            plt.show()

    def vnv_post(self, show=True):
        set_rcparams()
        c0, c1, c2 = get_default_color_palette()
        mk = get_default_markers()
        colors = [c0[0], c0[0], c0[1], c0[1], c0[2], c0[2], c0[3], c0[3], c0[4], c0[4]]
        markers= ['', mk[0], '', mk[1], '', mk[2], '', mk[3], '', mk[4]]
        ls = ['-', ':', '-', ':', '-', ':', '-', ':', '-', ':', '-', ':', '-', ':']
        fig, ax = plt.subplots(1, 2, figsize=(8.5, 3.5))

        LAB = [
          'Grass',
          'Grass (FD)',
          'MP$\&$M',
          'MP$\&$M (FD)',
          'MP$\&$M mod',
          'MP$\&$M mod (FD)',
          'VR',
          'VR (FD)',
          ]

        # Loop on cases:
        for i, case in enumerate(self.cases):
            RES = case['OUTPUT PARAMETERS']['RESULT FILE']
            file_ini = os.path.join(ROOT_DIR, 'RESU/'+RES+'ini.dat')
            file_fin = os.path.join(ROOT_DIR, 'RESU/'+RES+'fin.dat')
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
            fr= abs(u)/np.sqrt(9.81*h)

            if i==0:
                # plot h
                ax[0].plot(x, zb, color=colors[i], lw=0.75)
                #ax[0].fill_between(x, zf, -0.11, color='saddlebrown', alpha=0.65)
                #ax[0].fill_between(x, zb, zf, color='peru', alpha=0.45)
                #ax[0].fill_between(x, zb, zb+h, color='steelblue', alpha=0.25)
                ax[0].plot(x, zb+h, label=LAB[i], color=colors[i], 
                           marker=markers[i], markersize=2.5, markevery=5, lw=1., ls=ls[i])
                # plot b
                ax[1].plot(x, zb, label=LAB[i], color=colors[i], 
                           marker=markers[i], markersize=2.5, markevery=5, lw=1., ls=ls[i])
            else:
                # plot h
                ax[0].plot(x, zb, color=colors[i], lw=0.75)
                ax[0].plot(x, zb+h, label=LAB[i], color=colors[i], 
                           marker=markers[i], markersize=2.5, markevery=5, lw=1., ls=ls[i])
                # plot b
                ax[1].plot(x, zb, label=LAB[i], color=colors[i], 
                           marker=markers[i], markersize=2.5, markevery=5, lw=1., ls=ls[i])

        ax[0].legend()
        ax[0].set_ylabel("$z$ (m)")
        ax[0].set_xlabel("$x$ (m)")
        ax[0].grid()
        ax[0].set_xlim([3., 10.])
        ax[0].set_ylim([-0.01, 0.4])

        ax[1].legend()
        ax[1].set_ylabel("$z_b$ (m)")
        ax[1].set_xlabel("$x$ (m)")
        ax[1].grid()
        ax[1].set_xlim([3., 10.])
        ax[1].set_ylim([-0.011, 0.005])
        
        plt.savefig("FIG/bedload_dambreak_qs.png", dpi=300)
        if show:
            plt.show()

if __name__ == '__main__':

    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--run', default=False, action="store_true", help='run arrakis')
    parser.add_argument('-p', '--post', default=False, action="store_true", help='post processing')
    args = parser.parse_args()
    if args.run==False and args.post==False:
        no_args = True
    else:
        no_args = False

    # preparing vnv class
    vnv_case = VNV_bedload_dambreak_qs()
    vnv_case.vnv_pre()

    # Default run
    if no_args:
        vnv_case.vnv_run()
        vnv_case.vnv_post()

    # Custom run
    else:
        if args.run:
            vnv_case.vnv_run()
        if args.post:
            vnv_case.vnv_post()
