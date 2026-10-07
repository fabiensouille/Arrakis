#!/usr/bin/env python
import os
import numpy as np
import matplotlib.pylab as plt
from arrakis import *

ROOT_DIR = os.path.dirname(__file__)

class VNV_channel_constriction_10():

    def vnv_pre(self):
        self.name = "channel_constriction"
        self.cases = []

        # GLOBAL PARAMETERS
        Tf = 5400.
        Ks = 40.
        rhos = 2650.
        d50 = 2.e-3
        mpm_tausc = 0.053

        # analogy with MPM to get Ag
        g = 9.8065
        h0 = 0.066
        B0 = 0.064
        Rh = h0*B0/(2.*h0 + B0)
        R = (rhos/1000.)-1.
        grass_ag = 8.*np.sqrt(g)/(R*(Ks**3.)*np.sqrt(Rh))
        grass_mg = 3.
        print("Ag =", grass_ag)

        # Roe - Van Rijn
        case_inputs = read_yml_input_file(os.path.join(ROOT_DIR, "channel_constriction.yml"))
        case_inputs['TIME PARAMETERS']['DURATION'] = Tf
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'VR10'
        case_inputs['PHYSICAL PARAMETERS']['BOTTOM FRICTION COEFFICIENT'] = Ks
        case_inputs['NUMERICAL PARAMETERS']['ADVECTION SCHEME'] = 1
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['UNLOCK ADVANCED PARAMETERS'] = True
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['HYDROSTATIC RECONSTRUCTION'] = 4
        case_inputs['ADVANCED NUMERICAL PARAMETERS']['VARIABLE WIDTH SCHEME'] = 5
        case_inputs['SEDIMENT PARAMETERS']['NON ERODIBLE BED DAMPING REGION DEPTH'] = 0.01
        case_inputs['SEDIMENT PARAMETERS']['BEDLOAD FORMULA'] = 9
        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT DENSITY'] = rhos
        case_inputs['SEDIMENT PARAMETERS']['SEDIMENT D50'] = d50
        case_inputs['SEDIMENT PARAMETERS']['CRITICAL SHIELDS'] = mpm_tausc
        self.cases.append(case_inputs)


    def vnv_run(self, omp_num_threads):
        """
        Run 
        """
        rep = os.path.join(ROOT_DIR, "RESU")
        os.system("rm {}/*.dat".format(rep))
        for case in self.cases:
            run_arrakis(case, omp_num_threads=omp_num_threads)

    def vnv_check(self):
        """
        Check
        """
        for case in self.cases:
            RES = case['OUTPUT PARAMETERS']['RESULT FILE']
            file_fin = os.path.join(ROOT_DIR, 'RESU/'+RES+'fin.dat')
            file_ref = os.path.join(ROOT_DIR, 'REF/'+RES+'ref_{}.txt'.format(self.name))
            res = np.loadtxt(file_fin, delimiter=',', skiprows=2)
            ref = np.loadtxt(file_ref, delimiter=',', skiprows=2)
            np.testing.assert_allclose(res, ref, rtol=1e-7)

    def reset_ref(self):
        """
        Set reference file
        """
        for case in self.cases:
            RES = case['OUTPUT PARAMETERS']['RESULT FILE']
            file_fin = os.path.join(ROOT_DIR, 'RESU/'+RES+'fin.dat')
            file_ref = os.path.join(ROOT_DIR, 'REF/'+RES+'ref_{}.txt'.format(self.name))
            os.system("cp {} {}".format(file_fin, file_ref))

    def vnv_post(self, snapshot=1, show=True):
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
            file_ini = os.path.join(ROOT_DIR, 'RESU/'+RES+'ini.dat')
            file_fin = os.path.join(ROOT_DIR, 'RESU/'+RES+'_SNAP000{}.dat'.format(snapshot))
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
        
        # analytical solution
        if snapshot==1:
            ref_data = np.loadtxt('REF/ref10.txt', skiprows=2)
            ref_data_cr = np.loadtxt('REF/refCR_053_10.txt', skiprows=2)
        elif snapshot==2:
            ref_data = np.loadtxt('REF/ref30.txt', skiprows=2)
            ref_data_cr = np.loadtxt('REF/refCR_053_30.txt', skiprows=2)
        elif snapshot==3:
            ref_data = np.loadtxt('REF/ref60.txt', skiprows=2)
            ref_data_cr = np.loadtxt('REF/refCR_053_60.txt', skiprows=2)
        elif snapshot==4:
            ref_data = np.loadtxt('REF/ref90.txt', skiprows=2)
            ref_data_cr = np.loadtxt('REF/refCR_053_90.txt', skiprows=2)
        
        ax[1,0].plot(ref_data[:, 0], ref_data[:, 1], label="Obs.", color="k", ls='--', lw=1.)
        ax[1,0].plot(ref_data_cr[:, 0], ref_data_cr[:, 1], label="Capilla Roma $\&$ al.", color="r", ls='-.', lw=1.)
        
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
        ax[1,0].set_xlim([0.3 , 0.7])
        ax[1,0].set_ylim([0.04, 0.1])
        ax[1,0].grid()
        
        ax[2,1].legend(loc=2)
        ax[2,1].set_ylabel("$u$ (m/s)")
        ax[2,1].set_xlabel("$x$ (m)")
        ax[2,1].grid()
        ax[2,0].set_ylabel("$Fr$ (-)")
        ax[2,0].set_xlabel("$x$ (m)")
        ax[2,0].legend(loc=2)
        ax[2,0].grid()
        
        plt.savefig("FIG/{}_{}.png".format(self.name, snapshot), dpi=300)
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

    # preparing vnv class
    vnv_case = VNV_channel_constriction_10()
    vnv_case.vnv_pre()

    # Default run
    if no_args:
        vnv_case.vnv_run(args.ncsize)
        vnv_case.vnv_check()
        vnv_case.vnv_post(snapshot=1)
        vnv_case.vnv_post(snapshot=2)
        vnv_case.vnv_post(snapshot=3)
        vnv_case.vnv_post(snapshot=4)

    # Custom run
    else:
        if args.run:
            vnv_case.vnv_run(args.ncsize)
        if args.check:
            vnv_case.vnv_check()
        if args.post:
            vnv_case.vnv_post(snapshot=1)
            vnv_case.vnv_post(snapshot=2)
            vnv_case.vnv_post(snapshot=3)
            vnv_case.vnv_post(snapshot=4)
        if args.reset_ref:
            vnv_case.reset_ref()
