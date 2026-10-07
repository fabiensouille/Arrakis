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

        # case 1
        case_inputs = read_yml_input_file(os.path.join(ROOT_DIR, "dambreak.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'CASE-1'
        case_inputs['MESH PARAMETERS']['MESH FILE NAME'] = 'MESH_100'
        self.cases.append(case_inputs)

        # case 2 
        case_inputs = read_yml_input_file(os.path.join(ROOT_DIR, "dambreak.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'CASE-2'
        case_inputs['MESH PARAMETERS']['MESH FILE NAME'] = 'MESH_500'
        self.cases.append(case_inputs)
        
        # case 3
        case_inputs = read_yml_input_file(os.path.join(ROOT_DIR, "dambreak.yml"))
        case_inputs['OUTPUT PARAMETERS']['RESULT FILE'] = 'CASE-3'
        case_inputs['MESH PARAMETERS']['MESH FILE NAME'] = 'MESH_1000'
        self.cases.append(case_inputs)

    def vnv_run(self, omp_num_threads=1):
        """
        Run
        """
        rep = os.path.join(ROOT_DIR, "RESU")
        os.system("rm {}/*.dat".format(rep))
        # number of runs 
        nrepeat = 10
        # runs
        for case in self.cases:
            cputime = 0.
            for i in range(nrepeat):
                time_ini = time.time()
                run_arrakis(case, omp_num_threads=omp_num_threads)
                time_fin = time.time()
                cputime += time_fin-time_ini
            cputime /= nrepeat
            self.cputimes.append(cputime)

    def vnv_post(self, show=True):
        """
        Post
        """
        print("CPU time CASE-1 = ", self.cputimes[0])
        print("CPU time CASE-2 = ", self.cputimes[1])
        print("CPU time CASE-3 = ", self.cputimes[2])
        recputimes_seq = [0.015941262245178223, 0.060798025131225585, 0.19690101146697997]
        recputimes_par = [0.014632940292358398, 0.037631106376647946, 0.09895098209381104]

        # plot
        set_rcparams()
        c0, c1, c2 = get_default_color_palette()
        markers = get_default_markers()
        colors = c0 + c1 + c2
        fig, ax = plt.subplots(1, 1, figsize=(8.5, 8.))
        ax.plot([1, 2, 3], recputimes_seq, c='g', lw=0.5, ls='--', marker='s', label="SEQ")
        ax.plot([1, 2, 3], recputimes_par, c='r', lw=0.5, ls='--', marker='s', label="PAR (4)")
        ax.plot([1, 2, 3], self.cputimes, c='k', marker='o', label="PAR")
        ax.set_ylim([0., 0.2])
        plt.legend()
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
    vnv_case = VNV_benchmark()
    vnv_case.vnv_pre()

    # Default run
    if no_args:
        vnv_case.vnv_run(omp_num_threads=args.ncsize)
        vnv_case.vnv_post()

    # Custom run
    else:
        if args.run:
            vnv_case.vnv_run(omp_num_threads=args.ncsize)
        if args.post:
            vnv_case.vnv_post()
