#!/usr/bin/env python
import os
import numpy as np
import time
import matplotlib.pylab as plt
from arrakis import *

ROOT_DIR = os.path.dirname(__file__)

class VNV_benchmark():

    def vnv_pre(self):
        self.name = "vnv_benchmark"
        self.cases = []
        self.cputimes = []

        # HLL
        case_inputs = read_yml_input_file(os.path.join(ROOT_DIR, "channel_constriction.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'HLL-HR'
        case_inputs['NUMERICAL PARAMETERS']['ADVECTION SCHEME'] = 2
        case_inputs['NUMERICAL PARAMETERS']['HYDROSTATIC RECONSTRUCTION'] = 1
        case_inputs['NUMERICAL PARAMETERS']['VARIABLE WIDTH SCHEME'] = 2
        case_inputs['SEDIMENT PARAMETERS']['SVE WAVES APPROXIMATION'] = 1
        self.cases.append(case_inputs)

        # Roe 
        case_inputs = read_yml_input_file(os.path.join(ROOT_DIR, "channel_constriction.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'ROE-VC'
        case_inputs['NUMERICAL PARAMETERS']['ADVECTION SCHEME'] = 1
        case_inputs['NUMERICAL PARAMETERS']['HYDROSTATIC RECONSTRUCTION'] = 5
        case_inputs['NUMERICAL PARAMETERS']['VARIABLE WIDTH SCHEME'] = 5
        case_inputs['SEDIMENT PARAMETERS']['SVE WAVES APPROXIMATION'] = 0
        self.cases.append(case_inputs)

    def vnv_run(self):
        """
        Run 
        """
        rep = os.path.join(ROOT_DIR, "RESU")
        os.system("rm {}/*.dat".format(rep))
        # number of runs 
        nrepeat = 5
        # runs
        for case in self.cases:
            cputime = 0.
            for i in range(nrepeat):
                time_ini = time.time()
                run_arrakis(case)
                time_fin = time.time()
                cputime += time_fin-time_ini
            cputime /= nrepeat
            self.cputimes.append(cputime)

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

    def vnv_post(self, show=True):
        """
        Post
        """
        print("CPU time HLL = ", self.cputimes[0])
        print("CPU time ROE = ", self.cputimes[1])
        print("Ecart = ", self.cputimes[1]/self.cputimes[0])
        recputimes_v1p2 = [0.7342, 2.1497]
        recputimes_v1p3 = [0.2871, 0.7711]
        # plot
        set_rcparams()
        c0, c1, c2 = get_default_color_palette()
        markers = get_default_markers()
        colors = c0 + c1 + c2
        fig, ax = plt.subplots(1, 1, figsize=(8.5, 8.))
        ax.plot([1, 2], recputimes_v1p2, c='g', lw=0.5, ls='--', marker='s', label="v1p2")
        ax.plot([1, 2], recputimes_v1p3, c='r', lw=0.5, ls='--', marker='s', label="v1p3")
        ax.plot([1, 2], self.cputimes, c='k', marker='o')
        ax.set_ylim([0., 2.5])
        plt.legend()
        plt.show()

if __name__ == '__main__':

    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--run', default=False, action="store_true", help='run arrakis')
    parser.add_argument('-p', '--post', default=False, action="store_true", help='post processing')
    parser.add_argument('-c', '--check', default=False, action="store_true", help='check results')
    parser.add_argument('--reset-ref', default=False, action="store_true", help='WARNING: only use in extreme necessity')
    args = parser.parse_args()
    if args.run==False and args.post==False and args.check==False and args.reset_ref==False:
        no_args = True
    else:
        no_args = False

    # preparing vnv class
    vnv_case = VNV_benchmark()
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
